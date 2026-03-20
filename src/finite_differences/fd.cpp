/*
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2026 Benjamin F. Zwick
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors (alphabetically):
 *      George C. BOURANTAS
 *      Konstantinos A. MOUNTRIS
 *      Benjamin F. ZWICK
 */

#include "StreamVorti/finite_differences/fd.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace StreamVorti {

FiniteDiff::FiniteDiff(mfem::GridFunction &gf, int stencil_order)
    : gf_(&gf),
      stencil_order_(stencil_order)
{
    const mfem::FiniteElementSpace *fes = gf.FESpace();
    nnodes_ = fes->GetNDofs();
    dim_ = fes->GetMesh()->Dimension();

    if (stencil_order != 2 && stencil_order != 4)
    {
        MFEM_ABORT("FiniteDiff: stencil_order must be 2 or 4, got "
                   << stencil_order << ".");
    }
}

FiniteDiff::~FiniteDiff()
{}

void FiniteDiff::ExtractNodeCoordinates(
    std::vector<std::vector<double>> &coords) const
{
    const mfem::FiniteElementSpace *fes = gf_->FESpace();
    const mfem::Mesh *mesh = fes->GetMesh();

    coords.resize(dim_);
    for (int d = 0; d < dim_; ++d)
    {
        coords[d].resize(nnodes_);
    }

    // Use the GridFunction's own nodes via DofToVDof.  This is the correct
    // approach for any FE order because the GridFunction stores the
    // coordinates of every DOF (not just mesh vertices).
    //
    // For an H1 order-1 space, DOFs coincide with vertices but may have a
    // different ordering.  For higher-order spaces, there are interior DOFs
    // that have no corresponding mesh vertex.
    //
    // We need a GridFunction that holds the mesh node coordinates.  If the
    // mesh owns its own Nodes GridFunction (e.g., a higher-order mesh),
    // we use that.  Otherwise, we build one from the mesh vertices through
    // a temporary FE space.

    const mfem::GridFunction *mesh_nodes = mesh->GetNodes();

    if (mesh_nodes)
    {
        // Higher-order mesh: use the mesh's own node GridFunction
        const mfem::FiniteElementSpace *node_fes = mesh_nodes->FESpace();
        int node_ndofs = node_fes->GetNDofs();

        // The mesh nodes and the input gf may use different FE spaces.
        // If they have the same number of DOFs, map directly.
        if (node_ndofs == nnodes_)
        {
            for (int i = 0; i < nnodes_; ++i)
            {
                for (int d = 0; d < dim_; ++d)
                {
                    coords[d][i] = (*mesh_nodes)(node_fes->DofToVDof(i, d));
                }
            }
        }
        else
        {
            // Different FE space: project mesh nodes onto the input FE space
            // by evaluating vertex coordinates.  This handles order-1 input
            // FE space on a higher-order mesh.
            for (int i = 0; i < nnodes_; ++i)
            {
                // For order-1 H1 on any mesh, DOF i corresponds to vertex i
                if (i < mesh->GetNV())
                {
                    const double *vtx = mesh->GetVertex(i);
                    for (int d = 0; d < dim_; ++d)
                    {
                        coords[d][i] = vtx[d];
                    }
                }
                else
                {
                    MFEM_ABORT("FiniteDiff: DOF " << i
                               << " exceeds mesh vertex count "
                               << mesh->GetNV()
                               << ". Higher-order FD on higher-order meshes "
                               "requires matching FE spaces.");
                }
            }
        }
    }
    else
    {
        // Order-1 mesh without explicit nodes: DOFs are mesh vertices.
        if (nnodes_ != mesh->GetNV())
        {
            MFEM_ABORT("FiniteDiff: Number of DOFs (" << nnodes_
                       << ") does not match number of mesh vertices ("
                       << mesh->GetNV() << "). The FE space order may not "
                       "match the mesh order.");
        }
        for (int i = 0; i < nnodes_; ++i)
        {
            const double *vtx = mesh->GetVertex(i);
            for (int d = 0; d < dim_; ++d)
            {
                coords[d][i] = vtx[d];
            }
        }
    }
}

void FiniteDiff::FindUniqueCoordinates(const std::vector<double> &raw,
                                       std::vector<double> &unique,
                                       double tol_factor)
{
    unique = raw;
    std::sort(unique.begin(), unique.end());

    // Compute a mesh-aware absolute tolerance
    double span = unique.back() - unique.front();
    double tol = tol_factor * std::max(span, 1.0);

    unique.erase(std::unique(unique.begin(), unique.end(),
        [tol](double a, double b){ return std::abs(a - b) < tol; }),
        unique.end());
}

double FiniteDiff::VerifyUniformSpacing(const std::vector<double> &unique,
                                        const char *axis_name,
                                        double tol)
{
    int n = static_cast<int>(unique.size());
    if (n < 2)
    {
        MFEM_ABORT("FiniteDiff: Only " << n << " unique " << axis_name
                   << "-coordinates found. Need at least 2 for FD.");
    }

    double h_mean = (unique.back() - unique.front()) / (n - 1);

    double max_dev = 0.0;
    for (int i = 1; i < n; ++i)
    {
        double h_i = unique[i] - unique[i - 1];
        double dev = std::abs(h_i - h_mean) / h_mean;
        max_dev = std::max(max_dev, dev);
    }

    if (max_dev > tol)
    {
        MFEM_ABORT("FiniteDiff: Non-uniform grid spacing detected along "
                   << axis_name << "-axis. Maximum relative deviation: "
                   << max_dev << " (tolerance: " << tol << "). "
                   "Finite differences require a uniform structured grid. "
                   "Use DCPSE for non-uniform or unstructured meshes.");
    }

    return h_mean;
}

int FiniteDiff::CoordinateToIndex(double value,
                                  const std::vector<double> &unique,
                                  double tol)
{
    auto it = std::lower_bound(unique.begin(), unique.end(), value - tol);
    if (it != unique.end() && std::abs(*it - value) < tol)
    {
        return static_cast<int>(it - unique.begin());
    }
    return -1;
}

} // namespace StreamVorti
