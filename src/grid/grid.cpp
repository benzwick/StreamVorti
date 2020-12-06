/*
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
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
 */


#include "StreamVorti/grid/grid.hpp"


namespace StreamVorti {


Grid::Grid()
{}


Grid::Grid(const Grid &grid)
{
    // Assign by copying grid3D.
    *this = grid;
}


Grid::Grid(const mfem::GridFunction &gf)
{
    // See:
    // - https://github.com/mfem/mfem/issues/63#issuecomment-221646308
    // - https://github.com/mfem/mfem/issues/684#issuecomment-441872496

    mfem::Mesh *mesh = gf.FESpace()->GetMesh();
    const int dim = mesh->Dimension();
    // Note (number of nodes) > (number of vertices) when order > 1
    const int nNodes = gf.FESpace()->GetNDofs();
    const mfem::FiniteElementCollection *fec = gf.FESpace()->FEColl();

    // ASSERT gf is H1 (not L2, ND, RT, etc.)
    if (dynamic_cast<const mfem::H1_FECollection*>(fec) == nullptr)
    {
        MFEM_ABORT( "Grid function FE space is not H1." );
    }

    // ASSERT mesh dim == 1, 2 or 3
    if (dim < 1 || dim > 3)
    {
        MFEM_ABORT( "Mesh is " << dim << "D not 1D, 2D or 3D." );
    }

    // Create GridFunction with nodal coordinates
    mfem::FiniteElementSpace *fes = new mfem::FiniteElementSpace(mesh, fec, dim);
    mfem::GridFunction nodes(fes);
    mesh->GetNodes(nodes);

    double x, y, z;
    Node node;

    for (int i = 0; i < nNodes; ++i)
    {
        node.SetId(i);
        x = /*dim= 1*/ nodes(nodes.FESpace()->DofToVDof(i, 0));
        y = dim >= 2 ? nodes(nodes.FESpace()->DofToVDof(i, 1)) : 0.0;
        z = dim >= 3 ? nodes(nodes.FESpace()->DofToVDof(i, 2)) : 0.0;
        node.SetCoordinates(x, y, z);
        nodes_.emplace_back(node);
    }

    delete fes;
}


Grid::~Grid()
{}


void Grid::LoadFrom(const std::string &grid_filename)
{
    // Check if mesh filename is not empty.
    if (grid_filename.empty()) {
        throw std::invalid_argument(Logger::Error("Could not load grid. No grid filename was given.").c_str());
    }

    // Get the extension of the mesh filename.
    auto ext = grid_filename.substr(grid_filename.length()-4);

    // Clear the mesh containers.
    this->nodes_.clear();

    // Load the corresponding format.
    if (ext == ".inp") {
        AbaqusIO abaqus_io;
        abaqus_io.LoadMeshFrom(grid_filename.c_str());
        abaqus_io.LoadNodesIn(this->nodes_);
    }
    else {
        std::string error = Logger::Error("Could not load grid of unknown format. Expected [.inp] Check: ") + grid_filename;
        throw std::invalid_argument(error.c_str());
    }
}


bool Grid::operator == (const Grid &grid) const
{
    // Compare tetrahedral meshes for equality.
    return (this->nodes_ == grid.nodes_);
}


bool Grid::operator != (const Grid &grid) const
{
    // Compare grids for inequality.
    return !(*this == grid);
}


Grid & Grid::operator = (const Grid &grid)
{
    if (this != &grid) {
        // Assign nodes from grid.
        this->nodes_ = grid.nodes_;
    }

    return *this;
}



} //end of namespace StreamVorti
