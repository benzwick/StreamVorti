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
 */

/**
 * @file lisp_mesh.cpp
 * @brief Implementation of MFEM mesh bindings for Lisp
 */

#include "StreamVorti/lisp/lisp_mesh.hpp"
#include "StreamVorti/lisp/ecl_runtime.hpp"
#include "StreamVorti/lisp/lisp_bridge.hpp"

#include "mfem.hpp"

#include <fstream>
#include <cmath>
#include <set>

// ==================== C API Implementation ====================

extern "C" {

void* sv_make_cartesian_mesh_2d(int nx, int ny, int elem_type,
                                 double sx, double sy)
{
    mfem::Element::Type type;
    switch (elem_type) {
        case 2:
            type = mfem::Element::TRIANGLE;
            break;
        case 3:
        default:
            type = mfem::Element::QUADRILATERAL;
            break;
    }

    mfem::Mesh* mesh = new mfem::Mesh(nx, ny, type, false, sx, sy, false);
    return static_cast<void*>(mesh);
}

void* sv_make_cartesian_mesh_3d(int nx, int ny, int nz, int elem_type,
                                 double sx, double sy, double sz)
{
    mfem::Element::Type type;
    switch (elem_type) {
        case 4:
            type = mfem::Element::TETRAHEDRON;
            break;
        case 5:
        default:
            type = mfem::Element::HEXAHEDRON;
            break;
    }

    mfem::Mesh* mesh = new mfem::Mesh(nx, ny, nz, type, false, sx, sy, sz, false);
    return static_cast<void*>(mesh);
}

void* sv_load_mesh(const char* path, int format)
{
    try {
        // MFEM auto-detects format based on file extension
        // format parameter is for future use
        (void)format;

        mfem::Mesh* mesh = new mfem::Mesh(path, 1, 1);
        return static_cast<void*>(mesh);
    } catch (const std::exception& e) {
        return nullptr;
    }
}

int sv_save_mesh(void* mesh_ptr, const char* path, int format)
{
    if (!mesh_ptr) return -1;

    mfem::Mesh* mesh = static_cast<mfem::Mesh*>(mesh_ptr);

    try {
        std::ofstream ofs(path);
        if (!ofs) return -1;

        if (format == 2) {
            // VTK format
            mesh->PrintVTK(ofs);
        } else {
            // MFEM format (default)
            mesh->Print(ofs);
        }

        return 0;
    } catch (...) {
        return -1;
    }
}

void sv_free_mesh(void* mesh)
{
    delete static_cast<mfem::Mesh*>(mesh);
}

int sv_mesh_num_vertices(void* mesh)
{
    if (!mesh) return 0;
    return static_cast<mfem::Mesh*>(mesh)->GetNV();
}

int sv_mesh_num_elements(void* mesh)
{
    if (!mesh) return 0;
    return static_cast<mfem::Mesh*>(mesh)->GetNE();
}

int sv_mesh_num_boundary_elements(void* mesh)
{
    if (!mesh) return 0;
    return static_cast<mfem::Mesh*>(mesh)->GetNBE();
}

int sv_mesh_dimension(void* mesh)
{
    if (!mesh) return 0;
    return static_cast<mfem::Mesh*>(mesh)->Dimension();
}

int sv_mesh_space_dimension(void* mesh)
{
    if (!mesh) return 0;
    return static_cast<mfem::Mesh*>(mesh)->SpaceDimension();
}

void sv_mesh_get_vertices(void* mesh_ptr, double* coords, int* n)
{
    if (!mesh_ptr || !coords || !n) return;

    mfem::Mesh* mesh = static_cast<mfem::Mesh*>(mesh_ptr);
    int nv = mesh->GetNV();
    int dim = mesh->SpaceDimension();

    *n = nv * dim;

    for (int i = 0; i < nv; ++i) {
        const double* v = mesh->GetVertex(i);
        for (int d = 0; d < dim; ++d) {
            coords[i * dim + d] = v[d];
        }
    }
}

void sv_mesh_get_vertex(void* mesh_ptr, int vertex_id, double* coords)
{
    if (!mesh_ptr || !coords) return;

    mfem::Mesh* mesh = static_cast<mfem::Mesh*>(mesh_ptr);
    int dim = mesh->SpaceDimension();
    const double* v = mesh->GetVertex(vertex_id);

    for (int d = 0; d < dim; ++d) {
        coords[d] = v[d];
    }
}

void sv_mesh_get_element_vertices(void* mesh_ptr, int elem_id,
                                   int* vertices, int* n)
{
    if (!mesh_ptr || !vertices || !n) return;

    mfem::Mesh* mesh = static_cast<mfem::Mesh*>(mesh_ptr);
    mfem::Array<int> v;
    mesh->GetElementVertices(elem_id, v);

    *n = v.Size();
    for (int i = 0; i < v.Size(); ++i) {
        vertices[i] = v[i];
    }
}

int sv_mesh_get_element_type(void* mesh_ptr, int elem_id)
{
    if (!mesh_ptr) return -1;

    mfem::Mesh* mesh = static_cast<mfem::Mesh*>(mesh_ptr);
    return mesh->GetElementType(elem_id);
}

void sv_mesh_get_boundary_element_vertices(void* mesh_ptr, int bdr_elem_id,
                                            int* vertices, int* n)
{
    if (!mesh_ptr || !vertices || !n) return;

    mfem::Mesh* mesh = static_cast<mfem::Mesh*>(mesh_ptr);
    mfem::Array<int> v;
    mesh->GetBdrElementVertices(bdr_elem_id, v);

    *n = v.Size();
    for (int i = 0; i < v.Size(); ++i) {
        vertices[i] = v[i];
    }
}

int sv_mesh_get_boundary_attribute(void* mesh_ptr, int bdr_elem_id)
{
    if (!mesh_ptr) return -1;

    mfem::Mesh* mesh = static_cast<mfem::Mesh*>(mesh_ptr);
    return mesh->GetBdrAttribute(bdr_elem_id);
}

void sv_mesh_get_boundary_attributes(void* mesh_ptr, int* attrs, int* n)
{
    if (!mesh_ptr || !attrs || !n) return;

    mfem::Mesh* mesh = static_cast<mfem::Mesh*>(mesh_ptr);
    const mfem::Array<int>& bdr_attrs = mesh->bdr_attributes;

    *n = bdr_attrs.Size();
    for (int i = 0; i < bdr_attrs.Size(); ++i) {
        attrs[i] = bdr_attrs[i];
    }
}

int sv_mesh_set_boundary_attribute(void* mesh_ptr, int predicate_type,
                                    double value, double tolerance,
                                    int attribute)
{
    if (!mesh_ptr) return 0;

    mfem::Mesh* mesh = static_cast<mfem::Mesh*>(mesh_ptr);
    int count = 0;
    int dim = mesh->SpaceDimension();

    for (int i = 0; i < mesh->GetNBE(); ++i) {
        // Get center of boundary element
        mfem::Array<int> vertices;
        mesh->GetBdrElementVertices(i, vertices);

        double cx = 0.0, cy = 0.0, cz = 0.0;
        for (int v = 0; v < vertices.Size(); ++v) {
            const double* coords = mesh->GetVertex(vertices[v]);
            cx += coords[0];
            if (dim > 1) cy += coords[1];
            if (dim > 2) cz += coords[2];
        }
        cx /= vertices.Size();
        cy /= vertices.Size();
        cz /= vertices.Size();

        bool match = false;
        switch (predicate_type) {
            case 0: // x_eq
                match = std::abs(cx - value) < tolerance;
                break;
            case 1: // y_eq
                match = std::abs(cy - value) < tolerance;
                break;
            case 2: // z_eq
                match = std::abs(cz - value) < tolerance;
                break;
        }

        if (match) {
            mesh->SetBdrAttribute(i, attribute);
            ++count;
        }
    }

    return count;
}

void sv_mesh_refine(void* mesh_ptr, int times)
{
    if (!mesh_ptr) return;

    mfem::Mesh* mesh = static_cast<mfem::Mesh*>(mesh_ptr);
    for (int i = 0; i < times; ++i) {
        mesh->UniformRefinement();
    }
}

void* sv_make_grid_function(void* mesh_ptr, int order)
{
    if (!mesh_ptr) return nullptr;

    mfem::Mesh* mesh = static_cast<mfem::Mesh*>(mesh_ptr);
    int dim = mesh->Dimension();

    // Create FE collection and space (these are owned by GridFunction via raw pointers)
    mfem::H1_FECollection* fec = new mfem::H1_FECollection(order, dim);
    mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, fec, 1);
    mfem::GridFunction* gf = new mfem::GridFunction(fes);

    // Note: In production, we'd need to track ownership properly
    // For now, the GridFunction takes ownership of fes, but fec needs manual cleanup

    return static_cast<void*>(gf);
}

void sv_free_grid_function(void* gf)
{
    if (!gf) return;

    mfem::GridFunction* gfunc = static_cast<mfem::GridFunction*>(gf);
    // Get FESpace and FECollection before deleting
    mfem::FiniteElementSpace* fes =
        const_cast<mfem::FiniteElementSpace*>(gfunc->FESpace());
    const mfem::FiniteElementCollection* fec = fes->FEColl();

    delete gfunc;
    delete fes;
    delete fec;
}

int sv_grid_function_size(void* gf)
{
    if (!gf) return 0;
    return static_cast<mfem::GridFunction*>(gf)->Size();
}

void sv_grid_function_get_values(void* gf_ptr, double* values, int* n)
{
    if (!gf_ptr || !values || !n) return;

    mfem::GridFunction* gf = static_cast<mfem::GridFunction*>(gf_ptr);
    *n = gf->Size();

    for (int i = 0; i < *n; ++i) {
        values[i] = (*gf)(i);
    }
}

void sv_grid_function_set_values(void* gf_ptr, const double* values, int n)
{
    if (!gf_ptr || !values) return;

    mfem::GridFunction* gf = static_cast<mfem::GridFunction*>(gf_ptr);
    int size = std::min(n, gf->Size());

    for (int i = 0; i < size; ++i) {
        (*gf)(i) = values[i];
    }
}

void* sv_grid_function_get_fespace(void* gf)
{
    if (!gf) return nullptr;
    return const_cast<mfem::FiniteElementSpace*>(
        static_cast<mfem::GridFunction*>(gf)->FESpace());
}

void sv_fespace_get_dof_coords(void* fes_ptr, double* coords, int* n)
{
    if (!fes_ptr || !coords || !n) return;

    mfem::FiniteElementSpace* fes =
        static_cast<mfem::FiniteElementSpace*>(fes_ptr);
    mfem::Mesh* mesh = fes->GetMesh();

    int ndofs = fes->GetNDofs();
    int dim = mesh->SpaceDimension();

    *n = ndofs * dim;

    // Get DOF coordinates using GridFunction
    mfem::GridFunction nodes(fes);
    mesh->GetNodes(nodes);

    mfem::DenseMatrix dof_coords;
    fes->GetMesh()->GetNodes()->GetNodalValues(nodes, 1);

    // Extract coordinates for each DOF
    for (int i = 0; i < ndofs; ++i) {
        mfem::Vector pt;
        fes->GetMesh()->GetNode(i, pt);
        for (int d = 0; d < dim; ++d) {
            coords[i * dim + d] = pt(d);
        }
    }
}

} // extern "C"

// ==================== C++ API Implementation ====================

namespace StreamVorti {
namespace Lisp {

std::unique_ptr<mfem::Mesh>
MeshWrapper::makeCartesian2D(int nx, int ny, int elem_type,
                             double sx, double sy)
{
    mfem::Element::Type type;
    switch (elem_type) {
        case 2:
            type = mfem::Element::TRIANGLE;
            break;
        case 3:
        default:
            type = mfem::Element::QUADRILATERAL;
            break;
    }

    return std::make_unique<mfem::Mesh>(nx, ny, type, false, sx, sy, false);
}

std::unique_ptr<mfem::Mesh>
MeshWrapper::makeCartesian3D(int nx, int ny, int nz, int elem_type,
                             double sx, double sy, double sz)
{
    mfem::Element::Type type;
    switch (elem_type) {
        case 4:
            type = mfem::Element::TETRAHEDRON;
            break;
        case 5:
        default:
            type = mfem::Element::HEXAHEDRON;
            break;
    }

    return std::make_unique<mfem::Mesh>(nx, ny, nz, type, false, sx, sy, sz, false);
}

std::unique_ptr<mfem::Mesh>
MeshWrapper::loadFromFile(const std::string& path)
{
    return std::make_unique<mfem::Mesh>(path.c_str(), 1, 1);
}

void MeshWrapper::saveToFile(mfem::Mesh* mesh, const std::string& path,
                              int format)
{
    std::ofstream ofs(path);
    if (format == 2) {
        mesh->PrintVTK(ofs);
    } else {
        mesh->Print(ofs);
    }
}

std::vector<int>
MeshWrapper::findBoundaryElements(mfem::Mesh* mesh, BoundaryPredicate pred,
                                   double tolerance)
{
    std::vector<int> result;
    int dim = mesh->SpaceDimension();

    for (int i = 0; i < mesh->GetNBE(); ++i) {
        mfem::Array<int> vertices;
        mesh->GetBdrElementVertices(i, vertices);

        double cx = 0.0, cy = 0.0, cz = 0.0;
        for (int v = 0; v < vertices.Size(); ++v) {
            const double* coords = mesh->GetVertex(vertices[v]);
            cx += coords[0];
            if (dim > 1) cy += coords[1];
            if (dim > 2) cz += coords[2];
        }
        cx /= vertices.Size();
        cy /= vertices.Size();
        cz /= vertices.Size();

        if (pred(cx, cy, cz)) {
            result.push_back(i);
        }
    }

    return result;
}

int MeshWrapper::setBoundaryAttribute(mfem::Mesh* mesh,
                                       BoundaryPredicate pred,
                                       int attribute,
                                       double tolerance)
{
    auto elements = findBoundaryElements(mesh, pred, tolerance);
    for (int elem : elements) {
        mesh->SetBdrAttribute(elem, attribute);
    }
    return static_cast<int>(elements.size());
}

MeshWrapper::BoundaryPredicate
MeshWrapper::xEquals(double value, double tol)
{
    return [value, tol](double x, double, double) {
        return std::abs(x - value) < tol;
    };
}

MeshWrapper::BoundaryPredicate
MeshWrapper::yEquals(double value, double tol)
{
    return [value, tol](double, double y, double) {
        return std::abs(y - value) < tol;
    };
}

MeshWrapper::BoundaryPredicate
MeshWrapper::zEquals(double value, double tol)
{
    return [value, tol](double, double, double z) {
        return std::abs(z - value) < tol;
    };
}

void registerMeshFunctions()
{
    // Register mesh functions with ECL so they can be called from Lisp
    // This uses ECL's foreign function interface

    // The C functions are already exported and can be accessed via CFFI
    // This function is a placeholder for any additional registration needed

    // Example: Create Lisp wrapper functions
    Runtime::eval(R"(
        (defpackage :streamvorti.mesh
          (:use :cl)
          (:export #:make-cartesian-mesh-2d
                   #:make-cartesian-mesh-3d
                   #:load-mesh
                   #:save-mesh
                   #:mesh-num-vertices
                   #:mesh-num-elements
                   #:mesh-dimension))
    )");
}

} // namespace Lisp
} // namespace StreamVorti
