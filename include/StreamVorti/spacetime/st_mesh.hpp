/*
 * StreamVorti - Adaptive Space-Time Navier-Stokes Solver
 * Copyright (C) 2026 Benjamin F. Zwick
 *
 * Space-time mesh construction from spatial mesh.
 *
 * Constructs (d+1)-dimensional space-time slab meshes by extruding
 * a d-dimensional spatial mesh in the time direction.
 *
 * For 2D spatial problems: triangles -> tetrahedra or quads -> hexahedra
 * For 3D spatial problems: tetrahedra -> pentatopes or hexahedra -> tesseracts
 *
 * References:
 *   Karyofylli, Frings, Elgeti, Behr (2018). "Simplex space-time meshes
 *   in two-phase flow simulations."
 */

#ifndef STREAMVORTI_SPACETIME_MESH_HPP
#define STREAMVORTI_SPACETIME_MESH_HPP

#include "st_config.hpp"
#include "mfem.hpp"

namespace StreamVorti {
namespace SpaceTime {

/// Constructs a space-time slab mesh from a spatial mesh
///
/// Given a d-dimensional spatial mesh and a time interval [t_n, t_{n+1}],
/// creates a (d+1)-dimensional mesh where the extra dimension is time.
///
/// For simplex elements:
///   - 2D: triangles are split into tetrahedra
///   - 3D: tetrahedra are split into pentatopes (4-simplices)
///
/// For tensor-product elements:
///   - 2D: quads are extruded into hexahedra
///   - 3D: hexahedra are extruded into 4D hexahedra
class SpaceTimeMesh
{
public:
    /// Construct space-time slab mesh
    /// @param spatial_mesh  The d-dimensional spatial mesh
    /// @param t_start       Start time of the slab
    /// @param t_end         End time of the slab
    /// @param elem_type     Element type (Simplex or TensorProduct)
    /// @param n_time_steps  Number of time layers within the slab
    SpaceTimeMesh(mfem::Mesh &spatial_mesh,
                  double t_start, double t_end,
                  ElementType elem_type = ElementType::Simplex,
                  int n_time_steps = 1);

    ~SpaceTimeMesh();

    /// Get the space-time mesh
    mfem::Mesh *GetMesh() { return st_mesh_; }
    const mfem::Mesh *GetMesh() const { return st_mesh_; }

    /// Get the spatial mesh reference
    mfem::Mesh &GetSpatialMesh() { return spatial_mesh_; }

    /// Get the spatial dimension
    int GetSpatialDim() const { return spatial_dim_; }

    /// Get space-time dimension (spatial_dim + 1)
    int GetSpaceTimeDim() const { return spatial_dim_ + 1; }

    /// Get number of spatial vertices
    int GetNumSpatialVertices() const { return n_spatial_verts_; }

    /// Get number of time layers (including top and bottom)
    int GetNumTimeLayers() const { return n_time_layers_; }

    /// Get the time at a given layer
    double GetTimeAtLayer(int layer) const;

    /// Get boundary attribute for the bottom face (t = t_start)
    int GetBottomBdrAttr() const { return bottom_bdr_attr_; }

    /// Get boundary attribute for the top face (t = t_end)
    int GetTopBdrAttr() const { return top_bdr_attr_; }

    /// Get the spatial boundary attribute offset
    /// Spatial boundary attributes are shifted by this offset in ST mesh
    int GetSpatialBdrAttrOffset() const { return spatial_bdr_attr_offset_; }

    /// Map a spatial vertex index and time layer to a ST vertex index
    int SpatialToSTVertex(int spatial_vert, int time_layer) const
    {
        return time_layer * n_spatial_verts_ + spatial_vert;
    }

    /// Refine the space-time mesh uniformly
    void UniformRefinement();

    /// Refine selected elements
    void Refine(const mfem::Array<int> &elements_to_refine);

private:
    /// Build simplex space-time mesh from triangulated spatial mesh
    void BuildSimplexMesh(double t_start, double t_end, int n_time_steps);

    /// Build tensor-product space-time mesh from quad/hex spatial mesh
    void BuildTensorProductMesh(double t_start, double t_end, int n_time_steps);

    /// Set up boundary attributes for space-time mesh
    void SetupBoundaryAttributes();

    mfem::Mesh &spatial_mesh_;
    mfem::Mesh *st_mesh_;
    int spatial_dim_;
    int n_spatial_verts_;
    int n_time_layers_;
    ElementType elem_type_;

    int bottom_bdr_attr_;         ///< Attribute for t = t_n face
    int top_bdr_attr_;            ///< Attribute for t = t_{n+1} face
    int spatial_bdr_attr_offset_; ///< Offset for spatial boundary attrs

    std::vector<double> time_layers_; ///< Time values at each layer
};

#ifdef MFEM_USE_MPI
/// Parallel version of SpaceTimeMesh
class ParSpaceTimeMesh
{
public:
    ParSpaceTimeMesh(mfem::ParMesh &spatial_mesh,
                     double t_start, double t_end,
                     ElementType elem_type = ElementType::Simplex,
                     int n_time_steps = 1);

    ~ParSpaceTimeMesh();

    mfem::ParMesh *GetMesh() { return par_st_mesh_; }
    mfem::ParMesh &GetSpatialMesh() { return par_spatial_mesh_; }

    int GetBottomBdrAttr() const { return bottom_bdr_attr_; }
    int GetTopBdrAttr() const { return top_bdr_attr_; }
    int GetSpatialBdrAttrOffset() const { return spatial_bdr_attr_offset_; }

    void UniformRefinement();
    void Refine(const mfem::Array<int> &elements_to_refine);

private:
    mfem::ParMesh &par_spatial_mesh_;
    mfem::ParMesh *par_st_mesh_;
    int spatial_dim_;
    int bottom_bdr_attr_;
    int top_bdr_attr_;
    int spatial_bdr_attr_offset_;
};
#endif

} // namespace SpaceTime
} // namespace StreamVorti

#endif // STREAMVORTI_SPACETIME_MESH_HPP
