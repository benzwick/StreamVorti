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

#include "StreamVorti/finite_differences/par_fd.hpp"

#ifdef MFEM_USE_MPI

#include <iostream>

namespace StreamVorti {

ParFiniteDiff::ParFiniteDiff(mfem::ParGridFunction &gf, int stencil_order)
    : gf_(&gf),
      pfes_(gf.ParFESpace()),
      comm_(gf.ParFESpace()->GetComm()),
      rank_(0),
      nranks_(1),
      stencil_order_(stencil_order)
{
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &nranks_);

    if (stencil_order != 2 && stencil_order != 4)
    {
        MFEM_ABORT("ParFiniteDiff: stencil_order must be 2 or 4, got "
                   << stencil_order << ".");
    }
}

ParFiniteDiff::~ParFiniteDiff()
{}

} // namespace StreamVorti

#endif // MFEM_USE_MPI
