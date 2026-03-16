/*
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017 Konstantinos A. Mountris
 * Copyright (C) 2020-2025 Benjamin F. Zwick
 * Copyright (C) 2025 Weizheng Li
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
 *      Weizheng (Will) LI
 *      Konstantinos A. MOUNTRIS
 *      Benjamin F. ZWICK
 */

/*!
   \file streamvorti_par.hpp
   \brief Parallel StreamVorti header (MPI-enabled).

   Follows the same pattern as MFEM's pmesh.hpp: this file is only non-empty
   when MFEM_USE_MPI is defined. Include this header (instead of streamvorti.hpp)
   in translation units that use parallel features.
*/

#ifndef STREAMVORTI_PAR_HPP_
#define STREAMVORTI_PAR_HPP_

#include "StreamVorti/streamvorti.hpp"

#ifdef MFEM_USE_MPI

#include "StreamVorti/approximants/par_dcpse.hpp"
#include "StreamVorti/approximants/par_dcpse_2d.hpp"
#include "StreamVorti/support_domain/par_support_domain.hpp"

/**
 * @brief Initialize parallel DC PSE derivative operators
 *
 * Creates 2D or 3D parallel DC PSE object with ghost node exchange.
 * Performs neighbor search (including ghosts) and computes distributed
 * derivative matrices as HypreParMatrix. Timing information is printed to console.
 *
 * @param gf Parallel grid function containing nodal coordinates
 * @param dim Spatial dimension (2 or 3)
 * @param num_neighbors Number of neighbors for DC PSE stencil
 * @return Pointer to parallel DC PSE derivative object
 */
StreamVorti::ParDcpse2d* InitialiseParDCPSE(mfem::ParGridFunction& gf, int dim, int num_neighbors);

/**
 * @brief Extract velocity along a centerline for validation (parallel).
 *
 * Parallel version that correctly maps true DOF indices to local mesh vertices,
 * gathers data to rank 0 via MPI, and writes a single output file.
 *
 * @param u_velocity U-velocity field indexed by TRUE DOF (size = TrueVSize)
 * @param v_velocity V-velocity field indexed by TRUE DOF
 * @param pfes Parallel finite element space (provides DOF mapping + communicator)
 * @param filename Output file path (written by rank 0 only)
 * @param axis Fixed coordinate ('x' or 'y')
 * @param position Value of the fixed coordinate (e.g., 0.5 for x=0.5)
 * @param tol Tolerance for "on the line" detection
 */
void ExtractCenterline(const mfem::Vector& u_velocity,
                       const mfem::Vector& v_velocity,
                       mfem::ParFiniteElementSpace* pfes,
                       const std::string& filename,
                       char axis,
                       double position,
                       double tol = 0.01);

#endif // MFEM_USE_MPI
#endif // STREAMVORTI_PAR_HPP_
