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

#include "StreamVorti/approximants/par_dcpse.hpp"

#ifdef MFEM_USE_MPI

namespace StreamVorti {

ParDcpse::ParDcpse(mfem::ParGridFunction &gf, int num_neighbors)
    : ParSupportDomain(gf, num_neighbors),
      comm_(gf.ParFESpace()->GetComm()),
      rank_(0),
      nranks_(1)
{
    // Get MPI rank and size
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &nranks_);

    if (rank_ == 0) {
        std::cout << "ParDcpse: Initialized on " << nranks_ << " ranks" << std::endl;
    }
}

ParDcpse::~ParDcpse()
{
    // Derived classes handle cleanup of HypreParMatrix objects
}

} // namespace StreamVorti

#endif // MFEM_USE_MPI
