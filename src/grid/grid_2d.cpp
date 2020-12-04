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


#include "StreamVorti/grid/grid_2d.hpp"


namespace StreamVorti {


Grid2D::Grid2D()
{}


Grid2D::Grid2D(const Grid2D &grid2D)
{
    // Assign by copying grid3D.
    *this = grid2D;
}


Grid2D::Grid2D(const mfem::GridFunction &nodes)
// See
// - https://github.com/mfem/mfem/issues/63#issuecomment-221646308
// - mfem::GridFunction::ReorderByNodes()
{
	// ASSERT Ordering::byNODES (default)
	if (nodes.FESpace()->GetOrdering() != mfem::Ordering::byNODES)
	{
		MFEM_ABORT( "Node coordinates are not ordered by nodes." );
	}

	// ASSERT mesh dim == 2
	int dim = nodes.FESpace()->GetMesh()->Dimension();
	if (dim != 2)
	{
		MFEM_ABORT( "Mesh is " << dim << "D not 2D." );
	}

	const int nNodes = nodes.Size() / dim;
	double x, y;
	Node node;

	for (int i = 0; i < nNodes; ++i)
	{
		x = nodes(         i);
		y = nodes(nNodes + i);
		node.SetId(i);
		node.SetCoordinates(x, y, 0.);
		nodes_.emplace_back(node);
	}
}


Grid2D::~Grid2D()
{}


void Grid2D::LoadFrom(const std::string &grid_filename)
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


bool Grid2D::operator == (const Grid2D &grid2D) const
{
    // Compare tetrahedral meshes for equality.
    return (this->nodes_ == grid2D.nodes_);
}


bool Grid2D::operator != (const Grid2D &grid2D) const
{
    // Compare grids for inequality.
    return !(*this == grid2D);
}


Grid2D & Grid2D::operator = (const Grid2D &grid2D)
{
    if (this != &grid2D) {
        // Assign nodes from grid.
        this->nodes_ = grid2D.nodes_;
    }

    return *this;
}



} //end of namespace StreamVorti
