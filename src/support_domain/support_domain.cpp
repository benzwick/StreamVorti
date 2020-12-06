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
    typedef K::Point_3 Point_3d;
    typedef CGAL::Search_traits_3<K> TreeTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
    typedef Neighbor_search::Tree Tree;


    if (this->support_nodes_.empty()) {
        throw std::runtime_error(Logger::Error("Could not compute cutoff radiuses."
                                               " Should set the support nodes first").c_str());
    }

    this->cutoff_radiuses_.clear();
    this->cutoff_radiuses_.reserve(this->support_nodes_.size());

    std::vector<Point_3d> support_coords;
    for (auto node : this->support_nodes_) {
        support_coords.emplace_back(Point_3d(node.Coordinates().X(), node.Coordinates().Y(), node.Coordinates().Z()));
    }

    Tree tree(support_coords.begin(), support_coords.end());

    for (auto query_node : this->support_nodes_) {
        Point_3d query(query_node.Coordinates().X(), query_node.Coordinates().Y(), query_node.Coordinates().Z());

        // Initialize the search structure, and search all N points
        Neighbor_search search(tree, query, neighs_num);

        std::vector<double> distances;
        distances.reserve(neighs_num);

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
    typedef K::Point_3 Point_3d;
    typedef CGAL::Search_traits_3<K> TreeTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
    typedef Neighbor_search::Tree Tree;


    if (this->support_nodes_.empty()) {
        throw std::runtime_error(Logger::Error("Could not compute support radiuses."
                                               " Should set the support nodes first").c_str());
    }

    this->support_radiuses_.clear();
    this->support_radiuses_.reserve(this->support_nodes_.size());

    std::vector<Point_3d> support_coords;
    for (auto node : this->support_nodes_) {
        support_coords.emplace_back(Point_3d(node.Coordinates().X(), node.Coordinates().Y(), node.Coordinates().Z()));
    }

    Tree tree(support_coords.begin(), support_coords.end());

    for (auto query_node : this->support_nodes_) {
        Point_3d query(query_node.Coordinates().X(), query_node.Coordinates().Y(), query_node.Coordinates().Z());

        // Initialize the search structure, and search all N points
        Neighbor_search search(tree, query, neighs_num);

        std::vector<double> distances;
        distances.reserve(neighs_num);

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


const std::vector<std::vector<int> > SupportDomain::NeighborIndices()
{

    if (this->support_nodes_.empty()) {
        throw std::runtime_error(Logger::Error("Could not find neighbor nodes to given points."
                                               " Should set the support nodes first").c_str());
    }

    if (this->cutoff_radiuses_.empty()) {
        throw std::runtime_error(Logger::Error("Could not find neighbor nodes to given points."
                                               " Should compute the cutoff radiuses first").c_str());
    }


    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3d;
    typedef boost::tuple<Point_3d,int> Point_and_int;
    typedef CGAL::Search_traits_3<K> Traits_base;
    typedef CGAL::Search_traits_adapter<Point_and_int,
                                        CGAL::Nth_of_tuple_property_map<0, Point_and_int>, Traits_base> Traits;
    typedef CGAL::Kd_tree<Traits> Tree;
    typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;


    // Create container of neighbors indices for each node.
    std::vector<std::vector<int> > closest_nodes_ids;

    //Vectors to store coordinates and ids of point to be passed in tree tuples.
    std::vector<Point_3d> points;
    std::vector<int> points_indices;

    for (const auto &node : this->support_nodes_) {
        auto id = &node - &this->support_nodes_[0];
        points.emplace_back(Point_3d(node.Coordinates().X(), node.Coordinates().Y(), node.Coordinates().Z()));
        points_indices.emplace_back(id);
    }

    //Insert <point,id> tuples in the searching tree.
    Tree tree(boost::make_zip_iterator(boost::make_tuple(points.begin(), points_indices.begin() )),
              boost::make_zip_iterator(boost::make_tuple(points.end(), points_indices.end() )) );

    //Vector containing the domain nodes for a given grid node.
    std::vector<Point_and_int> domain_nodes;

    //Search domain_nodes for each grid node.
    for (auto &node : this->support_nodes_) {
        auto id = &node - &this->support_nodes_[0];

        Point_3d center(node.Coordinates().X(), node.Coordinates().Y(), node.Coordinates().Z());

        // Searching sphere.
        Fuzzy_sphere fs(center, this->cutoff_radiuses_[id]);

        //Neighbors search
        tree.search( std::back_inserter(domain_nodes), fs);

        // Vector to store domain_nodes indices.
        std::vector<int> neigh_indices;

        //Iterate over domain nodes.
        for (auto d_node : domain_nodes) {
            //Store domain nodes indices.
            neigh_indices.emplace_back(boost::get<1>(d_node));
        }

        closest_nodes_ids.emplace_back(neigh_indices);

        //Empty domain_nodes vector for next iteration.
        domain_nodes.clear();
    }


    // Return the indices of the closest nodes to each evaluation point.
    return closest_nodes_ids;

}


void SupportDomain::SaveNeighsToFile(const std::vector<std::vector<int> > &neighbor_ids,
                                     const std::string &filename) const
{
    //Initialize the path of the exporting file.
    std::string path = "";

    // Position of the last slash in the exporting file's name.
    std::size_t last_slash = filename.find_last_of("/\\");

    // Get the path directory of the exporting file name.
    if (last_slash != std::string::npos) {
        path = filename.substr(0, last_slash);
    }

    // Create the path's directory if it doesn't exist.
    boost::filesystem::path dir(path);
    if (!path.empty() && !boost::filesystem::exists(dir)) {
        boost::filesystem::create_directories(dir);
    }

    // Initialize the exporting file name's extension.
    std::string ext = "";

    // Search for the extension.
    if (filename.find_last_of(".") != std::string::npos) {
        ext = filename.substr(filename.find_last_of("."));
    }

    // Add .vtu extension before exporting if it's missing from the exporting file's name.
    std::string out_filename;
    if (ext != ".txt") { out_filename = filename + ".txt"; }
    else { out_filename = filename; }

    std::ofstream out(filename, std::ios::out | std::ios::trunc);

    for (auto &neighs : neighbor_ids) {
        for (auto &id : neighs) { out << id+1 << " "; }
        out << std::endl;
    }

    out.close();

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
