/*
 * ExplicitSim - Software for solving PDEs using explicit methods.
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
 *      Grand R. JOLDES
 *      Konstantinos A. MOUNTRIS
 */


#include "ExplicitSim/conditions/dirichlet.hpp"


namespace ExplicitSim {


Dirichlet::Dirichlet()
{}


Dirichlet::Dirichlet(const bool &x, const bool &y, const bool &z, const std::string &boundary_name)
{
    // Initialize the direction of the dirichlet condition.
    this->direction_ << 0, 0, 0;
    if (x) { this->direction_(0) = 1; }
    if (y) { this->direction_(1) = 1; }
    if (z) { this->direction_(2) = 1; }

    // Initialize the boundary id that the must have nodes where the dirichlet condition is applied.
    this->boundary_id_ = 0;

    // Initialize the boundary name that the must have nodes where the dirichlet condition is applied.
    this->boundary_name_ = boundary_name;
}


Dirichlet::~Dirichlet()
{}




} //end of namespace ExplicitSim
