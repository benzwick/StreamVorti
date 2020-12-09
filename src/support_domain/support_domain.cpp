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

SupportDomain::SupportDomain(const mfem::GridFunction &gf)
    :
    dim_(gf.FESpace()->GetMesh()->Dimension()),
    num_support_nodes_(gf.FESpace()->GetNDofs())
{
    mfem::Mesh *mesh = gf.FESpace()->GetMesh();
    const mfem::FiniteElementCollection *fec = gf.FESpace()->FEColl();

    // ASSERT gf is H1 (not L2, ND, RT, etc.)
    if (dynamic_cast<const mfem::H1_FECollection*>(fec) == nullptr)
    {
        MFEM_ABORT( "Grid function FE space is not H1." );
    }

    // ASSERT mesh 1 <= dim <= 3
    if (this->dim_ < 1 || this->dim_ > 3)
    {
        MFEM_ABORT( "Mesh is " << this->dim_ << "D not 1D, 2D or 3D." );
    }

    // Create GridFunctions with local and global nodal coordinates
    this->fespace_ = new mfem::FiniteElementSpace(mesh, fec, this->dim_);
    this->support_nodes_ = new mfem::GridFunction(this->fespace_);
    mesh->GetNodes(*this->support_nodes_);
    // TODO: this is needed for parallel implementation
    this->global_nodes_ = nullptr;
}


SupportDomain::SupportDomain(const mfem::ParGridFunction &gf, const mfem::Mesh &smesh)
    :
    dim_(gf.FESpace()->GetMesh()->Dimension()),
    num_support_nodes_(gf.FESpace()->GetTrueVSize())
{
    // TODO: these should be global not local support nodes
    // maybe pass the serial mesh to get the global nodes
    // maybe pass the global nodes as a GridFunction
    // maybe see `void ParGridFunction::SaveAsOne(std::ostream &out)`

    this->fespace_ = nullptr;
    this->support_nodes_ = nullptr;
    this->global_nodes_ = nullptr;
}


void SupportDomain::ComputeCutOffRadiuses(const std::size_t &neighs_num)
{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3d;
    typedef CGAL::Search_traits_3<K> TreeTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
    typedef Neighbor_search::Tree Tree;

    this->cutoff_radiuses_.clear();
    this->cutoff_radiuses_.reserve(this->num_support_nodes_);

    std::vector<Point_3d> support_coords;
    for (int i = 0; i < this->num_support_nodes_; ++i)
    {
        support_coords.emplace_back(this->SupportNodeAsPoint(i));
    }

    Tree tree(support_coords.begin(), support_coords.end());

    for (int i = 0; i < this->num_support_nodes_; ++i)
    {
        Point_3d query = this->SupportNodeAsPoint(i);

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

    this->support_radiuses_.clear();
    this->support_radiuses_.reserve(this->num_support_nodes_);

    std::vector<Point_3d> support_coords;
    for (int i = 0; i < this->num_support_nodes_; ++i)
    {
        support_coords.emplace_back(this->SupportNodeAsPoint(i));
    }

    Tree tree(support_coords.begin(), support_coords.end());

    for (int i = 0; i < this->num_support_nodes_; ++i)
    {
        Point_3d query = this->SupportNodeAsPoint(i);

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

    for (int id = 0; id < this->num_support_nodes_; ++id)
    {
        points.emplace_back(this->SupportNodeAsPoint(id));
        points_indices.emplace_back(id);
    }

    //Insert <point,id> tuples in the searching tree.
    Tree tree(boost::make_zip_iterator(boost::make_tuple(points.begin(), points_indices.begin() )),
              boost::make_zip_iterator(boost::make_tuple(points.end(), points_indices.end() )) );

    //Vector containing the domain nodes for a given grid node.
    std::vector<Point_and_int> domain_nodes;

    //Search domain_nodes for each grid node.
    std::vector<Point_3d> support_coords;
    for (int id = 0; id < this->num_support_nodes_; ++id)
    {
        Point_3d center = this->SupportNodeAsPoint(id);

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

} //end of namespace StreamVorti
