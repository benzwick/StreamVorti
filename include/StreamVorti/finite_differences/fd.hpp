/*
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017 Konstantinos A. Mountris
 * Copyright (C) 2020-2025 Benjamin F. Zwick
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

#ifndef STREAMVORTI_FINITE_DIFFERENCES_FD_HPP_
#define STREAMVORTI_FINITE_DIFFERENCES_FD_HPP_

#include "mfem.hpp"

#include "StreamVorti/derivative_operator.hpp"

#include <string>
#include <vector>

namespace StreamVorti {

/*!
 * \class FiniteDiff
 * \brief Finite difference derivatives base class for structured grids.
 *
 * Constructs sparse derivative matrices using standard finite difference
 * stencils on structured (Cartesian) MFEM meshes. Supports both 2nd-order
 * and 4th-order accurate stencils.
 *
 * For nodes at domain boundaries, one-sided stencils of matching order are
 * used to maintain uniform accuracy across the domain.
 *
 * The node coordinates are extracted from the GridFunction's
 * FiniteElementSpace, not from mesh vertices. This ensures correct behaviour
 * with higher-order (p>1) meshes used for visualization or isoparametric
 * elements.
 *
 * Grid validation:
 * - Verifies the mesh is structured (nx * ny [* nz] == nnodes)
 * - Verifies uniform spacing in each direction (max deviation < tol)
 * - Aborts with a diagnostic message if either check fails
 *
 * Usage:
 *
 *     GridFunction u_gf(fespace_h1);
 *     StreamVorti::FiniteDiff2d fd(u_gf);
 *     fd.Update();
 *     const SparseMatrix &dx = fd.ShapeFunctionDx();
 *     GridFunction dudx(fespace_h1);
 *     dx.Mult(u_gf, dudx);
 */
class FiniteDiff : public DerivativeOperator
{
public:
    /*!
     * \brief Constructor for structured MFEM grids.
     *
     * Extracts grid dimensions and spacing from the mesh associated with
     * the given GridFunction. The mesh must be a Cartesian structured grid.
     *
     * \param gf GridFunction whose FE space defines the node layout.
     * \param stencil_order Finite difference stencil order of accuracy
     *        (2 or 4). Default is 2.
     */
    FiniteDiff(mfem::GridFunction &gf, int stencil_order = 2);

    virtual ~FiniteDiff();

    virtual void Update() = 0;

    virtual void SaveDerivToFile(const std::string &deriv, const std::string &filename) const = 0;

    virtual const mfem::SparseMatrix & D(int i) const = 0;

    inline int NumNodes() const { return nnodes_; }

    inline int StencilOrder() const { return stencil_order_; }

protected:
    /*!
     * \brief Extract node coordinates from the GridFunction's FE space.
     *
     * Uses the GridFunction nodal values (not mesh vertices) so that
     * higher-order mesh representations are handled correctly.
     *
     * \param[out] coords Vector of coordinate vectors [x_coords, y_coords, ...]
     */
    void ExtractNodeCoordinates(std::vector<std::vector<double>> &coords) const;

    /*!
     * \brief Find unique sorted coordinate values along one axis.
     *
     * Uses a mesh-aware tolerance based on the coordinate range.
     *
     * \param raw    Raw coordinate values for all nodes
     * \param[out] unique Sorted unique values
     * \param tol_factor Relative tolerance factor (default 1e-10)
     */
    static void FindUniqueCoordinates(const std::vector<double> &raw,
                                      std::vector<double> &unique,
                                      double tol_factor = 1e-10);

    /*!
     * \brief Verify that grid spacing is uniform along one axis.
     *
     * Checks that max(|h_i - h_mean|) / h_mean < tol for all intervals.
     * Aborts with a diagnostic message if the grid is non-uniform.
     *
     * \param unique Sorted unique coordinate values along one axis
     * \param axis_name Name for diagnostic messages ("x", "y", or "z")
     * \param tol Relative tolerance for uniformity check (default 1e-6)
     * \return The uniform grid spacing h
     */
    static double VerifyUniformSpacing(const std::vector<double> &unique,
                                       const char *axis_name,
                                       double tol = 1e-6);

    /*!
     * \brief Map a raw coordinate to its grid index in the unique array.
     *
     * \param value  The coordinate value to look up
     * \param unique Sorted unique coordinate values
     * \param tol    Absolute tolerance for matching
     * \return Grid index, or -1 if not found
     */
    static int CoordinateToIndex(double value,
                                 const std::vector<double> &unique,
                                 double tol);

    mfem::GridFunction *gf_;
    int nnodes_;
    int dim_;
    int stencil_order_;
};

} // namespace StreamVorti

#endif // STREAMVORTI_FINITE_DIFFERENCES_FD_HPP_
