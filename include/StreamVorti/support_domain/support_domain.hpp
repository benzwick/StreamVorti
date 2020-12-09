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

/*!
   \file support_domain.hpp
   \brief SupportDomain class header file.
   \author Konstantinos A. Mountris
   \date 12/01/2018
*/

#ifndef STREAMVORTI_SUPPORT_DOMAIN_SUPPORT_DOMAIN_HPP_
#define STREAMVORTI_SUPPORT_DOMAIN_SUPPORT_DOMAIN_HPP_


#include "StreamVorti/elements/node.hpp"
#include "StreamVorti/vectors/vectors.hpp"
#include "StreamVorti/utilities/logger.hpp"

#include "mfem.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Fuzzy_sphere.h>
#include <boost/iterator/zip_iterator.hpp>

#include <boost/filesystem.hpp>

#include <iostream>
#include <fstream>
#include <list>
#include <cmath>
#include <string>
#include <stdexcept>
#include <exception>
#include <limits>


namespace StreamVorti {

/*!
 *  \addtogroup Approximants
 *  @{
 */


/*!
 * \class SupportDomain
 */

class SupportDomain
{
public:

    /*!
     * \brief Serial SupportDomain constructor.
     */
    SupportDomain(const mfem::GridFunction &gf);


    /*!
     * \brief Parallel SupportDomain constructor.
     */
    SupportDomain(const mfem::ParGridFunction &gf, const mfem::Mesh &smesh);


    /*!
     * \brief SupportDomain destructor.
     */
    virtual ~SupportDomain() {
        delete this->fespace_;
        delete this->support_nodes_;
        delete this->global_nodes_;
    };


    void ComputeCutOffRadiuses(const std::size_t &neighs_num);


    void ComputeSupportRadiuses(const std::size_t &neighs_num);


    inline const mfem::GridFunction & SupportNodes() const { return *this->support_nodes_; }


    inline CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 SupportNodeAsPoint(int id) {
        double x = 0.;
        double y = 0.;
        double z = 0.;
        if (dim_ > 0) x = (*support_nodes_)(fespace_->DofToVDof(id, 0));
        if (dim_ > 1) y = (*support_nodes_)(fespace_->DofToVDof(id, 1));
        if (dim_ > 2) z = (*support_nodes_)(fespace_->DofToVDof(id, 2));
        return CGAL::Exact_predicates_inexact_constructions_kernel::Point_3(x, y, z);
    }


    inline const std::vector<double> & CutoffRadiuses() const { return this->cutoff_radiuses_; }


    inline const std::vector<double> & SupportRadiuses() const { return this->support_radiuses_; }


    const std::vector<std::vector<int> > NeighborIndices();


    void SaveNeighsToFile(const std::vector<std::vector<int> > &neighbor_ids,
                          const std::string &filename) const;


private:
    /*! Dimension of the mesh */
    const int dim_;

    /*! Number of support nodes */
    const int num_support_nodes_;

    mfem::FiniteElementSpace *fespace_;

    /*! The local support nodes of the support domain. */
    mfem::GridFunction *support_nodes_;

    /*! The global support nodes of the support domain. */
    mfem::GridFunction *global_nodes_;

    std::vector<double> cutoff_radiuses_;

    /*! The support radiuses of the influence nodes of the support domain. */
    std::vector<double> support_radiuses_;
};


/*! @} End of Doxygen Groups*/
} //end of namespace StreamVorti

#endif //STREAMVORTI_SUPPORT_DOMAIN_SUPPORT_DOMAIN_HPP_
