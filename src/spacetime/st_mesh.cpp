/*
 * StreamVorti - Adaptive Space-Time Navier-Stokes Solver
 * Copyright (C) 2020-2026 Benjamin F. Zwick
 *
 * Space-time mesh construction implementation.
 */

#include <StreamVorti/spacetime/st_mesh.hpp>
#include <cassert>
#include <algorithm>
#include <iostream>

namespace StreamVorti {
namespace SpaceTime {

SpaceTimeMesh::SpaceTimeMesh(mfem::Mesh &spatial_mesh,
                             double t_start, double t_end,
                             ElementType elem_type,
                             int n_time_steps)
    : spatial_mesh_(spatial_mesh),
      st_mesh_(nullptr),
      spatial_dim_(spatial_mesh.Dimension()),
      n_spatial_verts_(spatial_mesh.GetNV()),
      n_time_layers_(n_time_steps + 1),
      elem_type_(elem_type),
      bottom_bdr_attr_(0),
      top_bdr_attr_(0),
      spatial_bdr_attr_offset_(0)
{
    assert(spatial_dim_ == 2 || spatial_dim_ == 3);
    assert(t_end > t_start);
    assert(n_time_steps >= 1);

    // Store time layer values
    time_layers_.resize(n_time_layers_);
    double dt = (t_end - t_start) / n_time_steps;
    for (int i = 0; i < n_time_layers_; ++i)
    {
        time_layers_[i] = t_start + i * dt;
    }

    if (elem_type == ElementType::Simplex)
    {
        BuildSimplexMesh(t_start, t_end, n_time_steps);
    }
    else
    {
        BuildTensorProductMesh(t_start, t_end, n_time_steps);
    }

    SetupBoundaryAttributes();
}

SpaceTimeMesh::~SpaceTimeMesh()
{
    delete st_mesh_;
}

double SpaceTimeMesh::GetTimeAtLayer(int layer) const
{
    assert(layer >= 0 && layer < n_time_layers_);
    return time_layers_[layer];
}

void SpaceTimeMesh::UniformRefinement()
{
    if (st_mesh_)
    {
        st_mesh_->UniformRefinement();
    }
}

void SpaceTimeMesh::Refine(const mfem::Array<int> &elements_to_refine)
{
    if (st_mesh_ && elements_to_refine.Size() > 0)
    {
        st_mesh_->GeneralRefinement(elements_to_refine);
    }
}

void SpaceTimeMesh::BuildSimplexMesh(double t_start, double t_end,
                                     int n_time_steps)
{
    // For a 2D spatial mesh, we create a 3D space-time mesh.
    // Each spatial triangle is extruded into a prism, then split into
    // 3 tetrahedra using Bey's algorithm (consistent diagonal choice).
    //
    // For a 3D spatial mesh, tetrahedra are extruded into pentahedra,
    // then split into pentatopes (4-simplices). However, MFEM only
    // supports up to 3D meshes, so for 3D+t we use a different approach
    // (multiple 3D spatial solves with DG coupling in time).

    const int st_dim = spatial_dim_ + 1;

    if (st_dim > 3)
    {
        // MFEM supports up to 3D. For 3D+t, we handle time-stepping
        // externally and only build 3D spatial meshes.
        std::cerr << "Warning: MFEM supports up to 3D meshes. "
                  << "For 3D spatial problems, time will be handled "
                  << "via semi-discrete DG coupling between time slabs."
                  << std::endl;
        // Just copy the spatial mesh
        st_mesh_ = new mfem::Mesh(spatial_mesh_, true);
        return;
    }

    // 2D spatial -> 3D space-time
    const int nv = n_spatial_verts_;
    const int ne = spatial_mesh_.GetNE();
    const int nbe = spatial_mesh_.GetNBE();
    const double dt = (t_end - t_start) / n_time_steps;

    // Total vertices: spatial vertices × time layers
    const int total_verts = nv * n_time_layers_;

    // For triangles -> tetrahedra: each tri-prism splits into 3 tets
    // For quads -> hexahedra split into tets (6 per hex) - but we're
    // in simplex mode, so spatial mesh should have triangles
    int tets_per_prism = 3;
    const int total_elems = ne * n_time_steps * tets_per_prism;

    // Boundary elements:
    // - Bottom face (t = t_start): ne spatial elements (triangles)
    // - Top face (t = t_end): ne spatial elements
    // - Lateral faces: nbe spatial boundary edges × n_time_steps,
    //   each extruded edge-rect splits into 2 triangles
    const int total_bdr = 2 * ne + nbe * n_time_steps * 2;

    st_mesh_ = new mfem::Mesh(st_dim, total_verts, total_elems, total_bdr);

    // Add vertices: for each time layer, copy spatial vertices with
    // additional time coordinate
    for (int tl = 0; tl < n_time_layers_; ++tl)
    {
        double t = time_layers_[tl];
        for (int iv = 0; iv < nv; ++iv)
        {
            const double *v = spatial_mesh_.GetVertex(iv);
            double st_coords[3] = {v[0], v[1], t};
            st_mesh_->AddVertex(st_coords);
        }
    }

    // Add elements: extrude each spatial element through each time step
    for (int ts = 0; ts < n_time_steps; ++ts)
    {
        for (int ie = 0; ie < ne; ++ie)
        {
            mfem::Array<int> verts;
            spatial_mesh_.GetElementVertices(ie, verts);
            int attr = spatial_mesh_.GetAttribute(ie);

            if (verts.Size() == 3)
            {
                // Triangle -> 3 tetrahedra (Bey's subdivision)
                // Bottom triangle: v0, v1, v2 at time layer ts
                // Top triangle: v0', v1', v2' at time layer ts+1
                int b0 = ts * nv + verts[0];
                int b1 = ts * nv + verts[1];
                int b2 = ts * nv + verts[2];
                int t0 = (ts + 1) * nv + verts[0];
                int t1 = (ts + 1) * nv + verts[1];
                int t2 = (ts + 1) * nv + verts[2];

                // Consistent diagonal: sort vertex indices to pick
                // a unique diagonal for the prism subdivision
                // Using the Dompierre et al. algorithm for consistent
                // simplex subdivision of prisms.

                // Tet 1: b0, b1, b2, t0
                st_mesh_->AddTet(b0, b1, b2, t0, attr);
                // Tet 2: b1, b2, t0, t1
                st_mesh_->AddTet(b1, b2, t0, t1, attr);
                // Tet 3: b2, t0, t1, t2
                st_mesh_->AddTet(b2, t0, t1, t2, attr);
            }
            else if (verts.Size() == 4)
            {
                // Quad -> 6 tetrahedra (standard hex-to-tet subdivision)
                int b0 = ts * nv + verts[0];
                int b1 = ts * nv + verts[1];
                int b2 = ts * nv + verts[2];
                int b3 = ts * nv + verts[3];
                int t0 = (ts + 1) * nv + verts[0];
                int t1 = (ts + 1) * nv + verts[1];
                int t2 = (ts + 1) * nv + verts[2];
                int t3 = (ts + 1) * nv + verts[3];

                // 5-tet subdivision of hexahedron (Kuhn's triangulation)
                st_mesh_->AddTet(b0, b1, b3, t0, attr);
                st_mesh_->AddTet(b1, b2, b3, t2, attr);
                st_mesh_->AddTet(b1, b3, t0, t2, attr);
                st_mesh_->AddTet(b1, t0, t1, t2, attr);
                st_mesh_->AddTet(b3, t0, t2, t3, attr);
            }
        }
    }

    // Add boundary elements

    // Bottom face (t = t_start): attribute for temporal inflow
    int max_spatial_bdr_attr = 0;
    for (int ibe = 0; ibe < nbe; ++ibe)
    {
        max_spatial_bdr_attr = std::max(max_spatial_bdr_attr,
                                        spatial_mesh_.GetBdrAttribute(ibe));
    }
    if (max_spatial_bdr_attr == 0) max_spatial_bdr_attr = 1;

    bottom_bdr_attr_ = max_spatial_bdr_attr + 1;
    top_bdr_attr_ = max_spatial_bdr_attr + 2;
    spatial_bdr_attr_offset_ = 0; // spatial attrs kept as-is

    // Bottom face: spatial elements at t = t_start
    for (int ie = 0; ie < ne; ++ie)
    {
        mfem::Array<int> verts;
        spatial_mesh_.GetElementVertices(ie, verts);

        if (verts.Size() == 3)
        {
            int b0 = verts[0]; // time layer 0
            int b1 = verts[1];
            int b2 = verts[2];
            st_mesh_->AddBdrTriangle(b0, b2, b1, bottom_bdr_attr_);
        }
        else if (verts.Size() == 4)
        {
            // Bottom quad -> 2 triangles
            int b0 = verts[0];
            int b1 = verts[1];
            int b2 = verts[2];
            int b3 = verts[3];
            st_mesh_->AddBdrTriangle(b0, b2, b1, bottom_bdr_attr_);
            st_mesh_->AddBdrTriangle(b0, b3, b2, bottom_bdr_attr_);
        }
    }

    // Top face: spatial elements at t = t_end
    int top_offset = n_time_steps * nv;
    for (int ie = 0; ie < ne; ++ie)
    {
        mfem::Array<int> verts;
        spatial_mesh_.GetElementVertices(ie, verts);

        if (verts.Size() == 3)
        {
            st_mesh_->AddBdrTriangle(top_offset + verts[0],
                                     top_offset + verts[1],
                                     top_offset + verts[2],
                                     top_bdr_attr_);
        }
        else if (verts.Size() == 4)
        {
            st_mesh_->AddBdrTriangle(top_offset + verts[0],
                                     top_offset + verts[1],
                                     top_offset + verts[2],
                                     top_bdr_attr_);
            st_mesh_->AddBdrTriangle(top_offset + verts[0],
                                     top_offset + verts[2],
                                     top_offset + verts[3],
                                     top_bdr_attr_);
        }
    }

    // Lateral faces: extrude spatial boundary edges
    for (int ts = 0; ts < n_time_steps; ++ts)
    {
        for (int ibe = 0; ibe < nbe; ++ibe)
        {
            mfem::Array<int> verts;
            spatial_mesh_.GetBdrElementVertices(ibe, verts);
            int bdr_attr = spatial_mesh_.GetBdrAttribute(ibe);

            if (verts.Size() == 2)
            {
                // Edge -> 2 triangles (quad face split)
                int b0 = ts * nv + verts[0];
                int b1 = ts * nv + verts[1];
                int t0 = (ts + 1) * nv + verts[0];
                int t1 = (ts + 1) * nv + verts[1];

                st_mesh_->AddBdrTriangle(b0, b1, t0, bdr_attr);
                st_mesh_->AddBdrTriangle(b1, t1, t0, bdr_attr);
            }
        }
    }

    st_mesh_->FinalizeTopology();
    st_mesh_->Finalize(true, true);
}

void SpaceTimeMesh::BuildTensorProductMesh(double t_start, double t_end,
                                           int n_time_steps)
{
    const int st_dim = spatial_dim_ + 1;

    if (st_dim > 3)
    {
        std::cerr << "Warning: MFEM supports up to 3D meshes. "
                  << "For 3D spatial problems with tensor-product elements, "
                  << "time will be handled via DG coupling between time slabs."
                  << std::endl;
        st_mesh_ = new mfem::Mesh(spatial_mesh_, true);
        return;
    }

    // 2D spatial -> 3D space-time with hexahedra
    const int nv = n_spatial_verts_;
    const int ne = spatial_mesh_.GetNE();
    const int nbe = spatial_mesh_.GetNBE();

    const int total_verts = nv * n_time_layers_;
    const int total_elems = ne * n_time_steps;
    const int total_bdr = 2 * ne + nbe * n_time_steps;

    st_mesh_ = new mfem::Mesh(st_dim, total_verts, total_elems, total_bdr);

    // Add vertices
    for (int tl = 0; tl < n_time_layers_; ++tl)
    {
        double t = time_layers_[tl];
        for (int iv = 0; iv < nv; ++iv)
        {
            const double *v = spatial_mesh_.GetVertex(iv);
            double st_coords[3] = {v[0], v[1], t};
            st_mesh_->AddVertex(st_coords);
        }
    }

    // Add elements: extrude quads into hexahedra
    for (int ts = 0; ts < n_time_steps; ++ts)
    {
        for (int ie = 0; ie < ne; ++ie)
        {
            mfem::Array<int> verts;
            spatial_mesh_.GetElementVertices(ie, verts);
            int attr = spatial_mesh_.GetAttribute(ie);

            if (verts.Size() == 4)
            {
                // Quad -> Hexahedron
                int b0 = ts * nv + verts[0];
                int b1 = ts * nv + verts[1];
                int b2 = ts * nv + verts[2];
                int b3 = ts * nv + verts[3];
                int t0 = (ts + 1) * nv + verts[0];
                int t1 = (ts + 1) * nv + verts[1];
                int t2 = (ts + 1) * nv + verts[2];
                int t3 = (ts + 1) * nv + verts[3];

                st_mesh_->AddHex(b0, b1, b2, b3, t0, t1, t2, t3, attr);
            }
            else if (verts.Size() == 3)
            {
                // Triangle -> Wedge/Prism (MFEM supports prisms)
                int b0 = ts * nv + verts[0];
                int b1 = ts * nv + verts[1];
                int b2 = ts * nv + verts[2];
                int t0 = (ts + 1) * nv + verts[0];
                int t1 = (ts + 1) * nv + verts[1];
                int t2 = (ts + 1) * nv + verts[2];

                st_mesh_->AddWedge(b0, b1, b2, t0, t1, t2, attr);
            }
        }
    }

    // Boundary elements
    int max_spatial_bdr_attr = 0;
    for (int ibe = 0; ibe < nbe; ++ibe)
    {
        max_spatial_bdr_attr = std::max(max_spatial_bdr_attr,
                                        spatial_mesh_.GetBdrAttribute(ibe));
    }
    if (max_spatial_bdr_attr == 0) max_spatial_bdr_attr = 1;

    bottom_bdr_attr_ = max_spatial_bdr_attr + 1;
    top_bdr_attr_ = max_spatial_bdr_attr + 2;
    spatial_bdr_attr_offset_ = 0;

    // Bottom face
    for (int ie = 0; ie < ne; ++ie)
    {
        mfem::Array<int> verts;
        spatial_mesh_.GetElementVertices(ie, verts);

        if (verts.Size() == 4)
        {
            st_mesh_->AddBdrQuad(verts[0], verts[3], verts[2], verts[1],
                                 bottom_bdr_attr_);
        }
        else if (verts.Size() == 3)
        {
            st_mesh_->AddBdrTriangle(verts[0], verts[2], verts[1],
                                     bottom_bdr_attr_);
        }
    }

    // Top face
    int top_offset = n_time_steps * nv;
    for (int ie = 0; ie < ne; ++ie)
    {
        mfem::Array<int> verts;
        spatial_mesh_.GetElementVertices(ie, verts);

        if (verts.Size() == 4)
        {
            st_mesh_->AddBdrQuad(top_offset + verts[0],
                                 top_offset + verts[1],
                                 top_offset + verts[2],
                                 top_offset + verts[3],
                                 top_bdr_attr_);
        }
        else if (verts.Size() == 3)
        {
            st_mesh_->AddBdrTriangle(top_offset + verts[0],
                                     top_offset + verts[1],
                                     top_offset + verts[2],
                                     top_bdr_attr_);
        }
    }

    // Lateral faces
    for (int ts = 0; ts < n_time_steps; ++ts)
    {
        for (int ibe = 0; ibe < nbe; ++ibe)
        {
            mfem::Array<int> verts;
            spatial_mesh_.GetBdrElementVertices(ibe, verts);
            int bdr_attr = spatial_mesh_.GetBdrAttribute(ibe);

            if (verts.Size() == 2)
            {
                // Edge -> Quad
                int b0 = ts * nv + verts[0];
                int b1 = ts * nv + verts[1];
                int t0 = (ts + 1) * nv + verts[0];
                int t1 = (ts + 1) * nv + verts[1];

                st_mesh_->AddBdrQuad(b0, b1, t1, t0, bdr_attr);
            }
        }
    }

    st_mesh_->FinalizeTopology();
    st_mesh_->Finalize(true, true);
}

void SpaceTimeMesh::SetupBoundaryAttributes()
{
    // Boundary attributes are already set during mesh construction.
    // Spatial boundary attributes are preserved as-is.
    // bottom_bdr_attr_ and top_bdr_attr_ are for temporal boundaries.
}

#ifdef MFEM_USE_MPI

ParSpaceTimeMesh::ParSpaceTimeMesh(mfem::ParMesh &spatial_mesh,
                                   double t_start, double t_end,
                                   ElementType elem_type,
                                   int n_time_steps)
    : par_spatial_mesh_(spatial_mesh),
      par_st_mesh_(nullptr),
      spatial_dim_(spatial_mesh.Dimension()),
      bottom_bdr_attr_(0),
      top_bdr_attr_(0),
      spatial_bdr_attr_offset_(0)
{
    // Build serial space-time mesh on this rank's partition
    mfem::Mesh local_spatial(spatial_mesh, false);
    SpaceTimeMesh serial_st(local_spatial, t_start, t_end,
                            elem_type, n_time_steps);

    bottom_bdr_attr_ = serial_st.GetBottomBdrAttr();
    top_bdr_attr_ = serial_st.GetTopBdrAttr();
    spatial_bdr_attr_offset_ = serial_st.GetSpatialBdrAttrOffset();

    // Create parallel mesh from the serial space-time mesh
    par_st_mesh_ = new mfem::ParMesh(spatial_mesh.GetComm(),
                                     *serial_st.GetMesh());
}

ParSpaceTimeMesh::~ParSpaceTimeMesh()
{
    delete par_st_mesh_;
}

void ParSpaceTimeMesh::UniformRefinement()
{
    if (par_st_mesh_)
    {
        par_st_mesh_->UniformRefinement();
    }
}

void ParSpaceTimeMesh::Refine(const mfem::Array<int> &elements_to_refine)
{
    if (par_st_mesh_ && elements_to_refine.Size() > 0)
    {
        par_st_mesh_->GeneralRefinement(elements_to_refine);
    }
}

#endif // MFEM_USE_MPI

} // namespace SpaceTime
} // namespace StreamVorti
