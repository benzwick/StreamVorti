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

#include "StreamVorti/approximants/dcpse.hpp"

namespace StreamVorti {

Dcpse::~Dcpse()
{}

Dcpse::Dcpse(mfem::GridFunction &gf,
             int CutoffRadAtNeighbor,
             int SupportRadAtNeighbor)
    : SupportDomain(gf)
{
    this->ComputeCutOffRadiuses(CutoffRadAtNeighbor);
    this->ComputeSupportRadiuses(SupportRadAtNeighbor);
}

} // namespace StreamVorti
