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
 * @file lisp_mesh.hpp
 * @brief MFEM mesh operations exposed to Lisp via CFFI
 *
 * Provides C-style functions for mesh creation, manipulation, and querying
 * that can be called from Common Lisp through CFFI.
 */

#ifndef STREAMVORTI_LISP_LISP_MESH_HPP_
#define STREAMVORTI_LISP_LISP_MESH_HPP_

#include <cstdint>

// C-style API for CFFI
#ifdef __cplusplus
extern "C" {
#endif

// ==================== Mesh Creation ====================

/**
 * @brief Create a 2D Cartesian mesh
 *
 * @param nx Number of elements in x direction
 * @param ny Number of elements in y direction
 * @param elem_type Element type: 2=Triangle, 3=Quadrilateral
 * @param sx Size in x direction (default 1.0)
 * @param sy Size in y direction (default 1.0)
 * @return Opaque pointer to MFEM Mesh
 */
void* sv_make_cartesian_mesh_2d(int nx, int ny, int elem_type,
                                 double sx, double sy);

/**
 * @brief Create a 3D Cartesian mesh
 *
 * @param nx Number of elements in x direction
 * @param ny Number of elements in y direction
 * @param nz Number of elements in z direction
 * @param elem_type Element type: 4=Tetrahedron, 5=Hexahedron
 * @param sx Size in x direction
 * @param sy Size in y direction
 * @param sz Size in z direction
 * @return Opaque pointer to MFEM Mesh
 */
void* sv_make_cartesian_mesh_3d(int nx, int ny, int nz, int elem_type,
                                 double sx, double sy, double sz);

/**
 * @brief Load mesh from file
 *
 * @param path Path to mesh file
 * @param format Format code: 0=auto, 1=mfem, 2=gmsh, 3=vtk, 4=netgen
 * @return Opaque pointer to MFEM Mesh, or NULL on error
 */
void* sv_load_mesh(const char* path, int format);

/**
 * @brief Save mesh to file
 *
 * @param mesh Pointer to mesh
 * @param path Output file path
 * @param format Format code: 1=mfem, 2=vtk
 * @return 0 on success, -1 on error
 */
int sv_save_mesh(void* mesh, const char* path, int format);

/**
 * @brief Free mesh memory
 */
void sv_free_mesh(void* mesh);

// ==================== Mesh Queries ====================

/**
 * @brief Get number of vertices/nodes
 */
int sv_mesh_num_vertices(void* mesh);

/**
 * @brief Get number of elements
 */
int sv_mesh_num_elements(void* mesh);

/**
 * @brief Get number of boundary elements
 */
int sv_mesh_num_boundary_elements(void* mesh);

/**
 * @brief Get mesh dimension (2 or 3)
 */
int sv_mesh_dimension(void* mesh);

/**
 * @brief Get spatial dimension (embedding dimension)
 */
int sv_mesh_space_dimension(void* mesh);

/**
 * @brief Get vertex coordinates
 *
 * @param mesh Pointer to mesh
 * @param coords Output array (must be pre-allocated: num_vertices * space_dim)
 * @param n Output: number of coordinates written
 */
void sv_mesh_get_vertices(void* mesh, double* coords, int* n);

/**
 * @brief Get single vertex coordinates
 *
 * @param mesh Pointer to mesh
 * @param vertex_id Vertex index (0-based)
 * @param coords Output array (must be pre-allocated: space_dim)
 */
void sv_mesh_get_vertex(void* mesh, int vertex_id, double* coords);

/**
 * @brief Get element vertices
 *
 * @param mesh Pointer to mesh
 * @param elem_id Element index (0-based)
 * @param vertices Output array for vertex indices
 * @param n Output: number of vertices
 */
void sv_mesh_get_element_vertices(void* mesh, int elem_id,
                                   int* vertices, int* n);

/**
 * @brief Get element type
 *
 * @param mesh Pointer to mesh
 * @param elem_id Element index (0-based)
 * @return Element type code
 */
int sv_mesh_get_element_type(void* mesh, int elem_id);

/**
 * @brief Get boundary element vertices
 */
void sv_mesh_get_boundary_element_vertices(void* mesh, int bdr_elem_id,
                                            int* vertices, int* n);

/**
 * @brief Get boundary attribute
 */
int sv_mesh_get_boundary_attribute(void* mesh, int bdr_elem_id);

/**
 * @brief Get all boundary attributes (unique values)
 *
 * @param mesh Pointer to mesh
 * @param attrs Output array
 * @param n Output: number of attributes
 */
void sv_mesh_get_boundary_attributes(void* mesh, int* attrs, int* n);

// ==================== Mesh Modification ====================

/**
 * @brief Set boundary attribute for elements matching predicate
 *
 * @param mesh Pointer to mesh
 * @param predicate_type 0=x_eq, 1=y_eq, 2=z_eq, 3=custom
 * @param value Value for predicate (coordinate value for eq predicates)
 * @param tolerance Tolerance for floating point comparison
 * @param attribute Attribute value to set
 * @return Number of boundary elements modified
 */
int sv_mesh_set_boundary_attribute(void* mesh, int predicate_type,
                                    double value, double tolerance,
                                    int attribute);

/**
 * @brief Refine mesh uniformly
 *
 * @param mesh Pointer to mesh
 * @param times Number of refinement iterations
 */
void sv_mesh_refine(void* mesh, int times);

// ==================== GridFunction Creation ====================

/**
 * @brief Create a GridFunction on the mesh
 *
 * @param mesh Pointer to mesh
 * @param order Polynomial order (1 = linear)
 * @return Opaque pointer to GridFunction
 */
void* sv_make_grid_function(void* mesh, int order);

/**
 * @brief Free GridFunction memory
 */
void sv_free_grid_function(void* gf);

/**
 * @brief Get number of DOFs
 */
int sv_grid_function_size(void* gf);

/**
 * @brief Get GridFunction values
 */
void sv_grid_function_get_values(void* gf, double* values, int* n);

/**
 * @brief Set GridFunction values
 */
void sv_grid_function_set_values(void* gf, const double* values, int n);

// ==================== FiniteElementSpace ====================

/**
 * @brief Get FiniteElementSpace from GridFunction
 */
void* sv_grid_function_get_fespace(void* gf);

/**
 * @brief Get DOF coordinates from FiniteElementSpace
 *
 * @param fes Pointer to FiniteElementSpace
 * @param coords Output array (num_dofs * dim)
 * @param n Output: number of coordinates
 */
void sv_fespace_get_dof_coords(void* fes, double* coords, int* n);

#ifdef __cplusplus
}
#endif

// C++ API
#ifdef __cplusplus

#include <memory>
#include <string>
#include <vector>
#include <functional>

namespace mfem {
    class Mesh;
    class GridFunction;
    class FiniteElementSpace;
}

namespace StreamVorti {
namespace Lisp {

/**
 * @class MeshWrapper
 * @brief C++ wrapper for MFEM Mesh with Lisp integration
 */
class MeshWrapper {
public:
    /**
     * @brief Create 2D Cartesian mesh
     */
    static std::unique_ptr<mfem::Mesh>
    makeCartesian2D(int nx, int ny, int elem_type = 3,
                    double sx = 1.0, double sy = 1.0);

    /**
     * @brief Create 3D Cartesian mesh
     */
    static std::unique_ptr<mfem::Mesh>
    makeCartesian3D(int nx, int ny, int nz, int elem_type = 5,
                    double sx = 1.0, double sy = 1.0, double sz = 1.0);

    /**
     * @brief Load mesh from file
     */
    static std::unique_ptr<mfem::Mesh>
    loadFromFile(const std::string& path);

    /**
     * @brief Save mesh to file
     */
    static void saveToFile(mfem::Mesh* mesh, const std::string& path,
                           int format = 1);

    /**
     * @brief Find boundary elements matching predicate
     *
     * Predicate receives (x, y, z) coordinates of boundary element center
     */
    using BoundaryPredicate = std::function<bool(double, double, double)>;

    static std::vector<int>
    findBoundaryElements(mfem::Mesh* mesh, BoundaryPredicate pred,
                         double tolerance = 1e-10);

    /**
     * @brief Set boundary attribute for matching elements
     */
    static int setBoundaryAttribute(mfem::Mesh* mesh,
                                    BoundaryPredicate pred,
                                    int attribute,
                                    double tolerance = 1e-10);

    /**
     * @brief Create predicate for x = value
     */
    static BoundaryPredicate xEquals(double value, double tol = 1e-10);

    /**
     * @brief Create predicate for y = value
     */
    static BoundaryPredicate yEquals(double value, double tol = 1e-10);

    /**
     * @brief Create predicate for z = value
     */
    static BoundaryPredicate zEquals(double value, double tol = 1e-10);
};

/**
 * @brief Register mesh functions with ECL
 *
 * Call after ECL initialization to make mesh functions available from Lisp.
 */
void registerMeshFunctions();

} // namespace Lisp
} // namespace StreamVorti

#endif // __cplusplus

#endif // STREAMVORTI_LISP_LISP_MESH_HPP_
