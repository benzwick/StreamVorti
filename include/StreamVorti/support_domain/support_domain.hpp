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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
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
     * \brief SupportDomain constructor.
     */
    SupportDomain();


    /*!
     * \brief SupportDomain destructor.
     */
    virtual ~SupportDomain();


    inline void SetSupportNodes(const std::vector<Node> &support_nodes) { this->support_nodes_ = support_nodes; }


    void ComputeCutOffRadiuses(const std::size_t &neighs_num);


    void ComputeSupportRadiuses(const std::size_t &neighs_num);


    inline const std::vector<Node> & SupportNodes() const { return this->support_nodes_; }


    inline const std::vector<double> & CutoffRadiuses() const { return this->cutoff_radiuses_; }


    inline const std::vector<double> & SupportRadiuses() const { return this->support_radiuses_; }


    const std::vector<std::vector<int> > NeighborIndices();


    void SaveNeighsToFile(const std::vector<std::vector<int> > &neighbor_ids,
                          const std::string &filename) const;


    int MinSupportNodesIn(const std::vector<std::vector<int> > &neighbor_ids) const;


    int MaxSupportNodesIn(const std::vector<std::vector<int> > &neighbor_ids) const;


private:
    std::vector<Node> support_nodes_;                /*!< The support nodes of the support domain. */

    std::vector<double> cutoff_radiuses_;

    std::vector<double> support_radiuses_;           /*!< The support radiuses of the influence nodes of the support domain. */
};


/*! @} End of Doxygen Groups*/
} //end of namespace StreamVorti

#endif //STREAMVORTI_SUPPORT_DOMAIN_SUPPORT_DOMAIN_HPP_
