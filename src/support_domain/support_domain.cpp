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

#include "StreamVorti/support_domain/support_domain.hpp"


namespace StreamVorti {

SupportDomain::SupportDomain()
{}


SupportDomain::~SupportDomain()
{}


void SupportDomain::ComputeCutOffRadiuses(const std::size_t &neighs_num)
{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_2 Point_2d;
    typedef CGAL::Search_traits_2<K> TreeTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
    typedef Neighbor_search::Tree Tree;


    if (this->support_nodes_.empty()) {
        throw std::runtime_error(Logger::Error("Could not compute cutoff radiuses."
                                               " Should set the support nodes first").c_str());
    }

    this->cutoff_radiuses_.clear();
    this->cutoff_radiuses_.reserve(this->support_nodes_.size());

    std::vector<Point_2d> support_coords;
    for (auto node : this->support_nodes_) {
        support_coords.emplace_back(Point_2d(node.Coordinates().X(), node.Coordinates().Y()));
    }

    Tree tree(support_coords.begin(), support_coords.end());

    for (auto query_node : this->support_nodes_) {
        Point_2d query(query_node.Coordinates().X(), query_node.Coordinates().Y());

        // Initialize the search structure, and search all N points
        Neighbor_search search(tree, query, neighs_num);

        std::vector<double> distances;

        // report the N nearest neighbors and their distance
        // This should sort all N points by increasing distance from origin
        for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
          distances.emplace_back(std::sqrt(it->second));
        }

        // Sort distances.
        std::sort(distances.begin(), distances.end());

        this->cutoff_radiuses_.emplace_back(distances.back());
    }


}


void SupportDomain::ComputeSupportRadiuses(const std::size_t &neighs_num)
{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_2 Point_2d;
    typedef CGAL::Search_traits_2<K> TreeTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
    typedef Neighbor_search::Tree Tree;


    if (this->support_nodes_.empty()) {
        throw std::runtime_error(Logger::Error("Could not compute support radiuses."
                                               " Should set the support nodes first").c_str());
    }

    this->support_radiuses_.clear();
    this->support_radiuses_.reserve(this->support_nodes_.size());

    std::vector<Point_2d> support_coords;
    for (auto node : this->support_nodes_) {
        support_coords.emplace_back(Point_2d(node.Coordinates().X(), node.Coordinates().Y()));
    }

    Tree tree(support_coords.begin(), support_coords.end());

    for (auto query_node : this->support_nodes_) {
        Point_2d query(query_node.Coordinates().X(), query_node.Coordinates().Y());

        // Initialize the search structure, and search all N points
        Neighbor_search search(tree, query, neighs_num);

        std::vector<double> distances;

        // report the N nearest neighbors and their distance
        // This should sort all N points by increasing distance from origin
        for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
          distances.emplace_back(std::sqrt(it->second));
        }

        // Sort distances.
        std::sort(distances.begin(), distances.end());

        this->support_radiuses_.emplace_back(distances.back());
    }


}


const std::vector<std::vector<int> > SupportDomain::CgalNeighborIdsTo(const std::vector<Vec3<double> > &eval_nodes_coords,
                                                                      const size_t &neigs_num)
{


//    if (this->support_nodes_.empty()) {
//        throw std::runtime_error(Logger::Error("Could not find neighbor nodes to given points."
//                                               " Should set the support nodes first").c_str());
//    }



//    // Return the indices of the closest nodes to each evaluation point.
//    return closest_nodes_ids;

}


int SupportDomain::MinSupportNodesIn(const std::vector<std::vector<int> > &neighbor_ids) const
{
    // Initialize minimum number of support nodes.
    int min_support_nodes_ = std::numeric_limits<int>::max();

    // Compute the minimum number of support nodes.
    for (auto &neighs : neighbor_ids) {
        if (static_cast<int>(neighs.size()) < min_support_nodes_) {
            min_support_nodes_ = static_cast<int>(neighs.size());
        }
    }

    return min_support_nodes_;
}


int SupportDomain::MaxSupportNodesIn(const std::vector<std::vector<int> > &neighbor_ids) const
{
    // Initialize maximum number of support nodes.
    int max_support_nodes_ = 0;

    // Compute the maximum number of support nodes.
    for (auto &neighs : neighbor_ids) {
        if (static_cast<int>(neighs.size()) > max_support_nodes_) {
            max_support_nodes_ = static_cast<int>(neighs.size());
        }
    }

    return max_support_nodes_;
}


} //end of namespace StreamVorti
