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


#include "ExplicitSim/conditions/ebciem.hpp"


namespace ExplicitSim {


Ebciem::Ebciem()
{}


Ebciem::~Ebciem()
{}


void Ebciem::SetNodalMmlsProperties(const std::string &basis_func_type, const bool &use_exact_derivatives)
{
    this->nodal_mmls_.SetBasisFunctionType(basis_func_type);
    this->nodal_mmls_.SetExactDerivativesMode(use_exact_derivatives);
}


void Ebciem::ComputeCorrectionMatrices(const WeakModel3D &model, const SupportDomain &support,
                                       const std::vector<Dirichlet> &dirichlet_conds,
                                       const std::vector<Loading> &loading_conds, bool is_simplified)
{
    // Check if dirichlet conditions are initialized.
    if (dirichlet_conds.empty()) {
        throw std::invalid_argument(Logger::Error("Could not apply EBCIEM correction. "
                                                  "The Dirichlet conditions container is empty.").c_str());
    }

    // Check if loading conditions are initialized.
    if (loading_conds.empty()) {
        throw std::invalid_argument(Logger::Error("Could not apply EBCIEM correction. "
                                                  "The Loading conditions container is empty.").c_str());
    }


    // Find influence nodes indices of the model's mesh nodes.
    auto neighs_ids = support.CgalClosestNodesIdsTo(model.TetrahedralMesh().NodeCoordinates());

    // Compute shape functions and derivatives.
    this->nodal_mmls_.ComputeShFuncAndDerivs(model.TetrahedralMesh().Nodes(), model.TetrahedralMesh().NodeCoordinates(),
                                             neighs_ids, support.InfluenceNodesRadiuses());

    // Initialize inverse mass matrix triplet.
    std::vector<Eigen::Triplet<double> > inv_mass_trip;
    inv_mass_trip.reserve(model.TetrahedralMesh().NodesNum());

    // Get inverse matrix values for each node.
    for (const auto &nodal_mass : model.Mass()) {
        auto node_id = &nodal_mass - &model.Mass()[0];
        inv_mass_trip.emplace_back(Eigen::Triplet<double>(node_id, node_id, 1./nodal_mass));
    }

    // Create inverse mass matrix.
    Eigen::SparseMatrix<double> inv_mass(model.TetrahedralMesh().NodesNum(), model.TetrahedralMesh().NodesNum());
    inv_mass.setFromTriplets(inv_mass_trip.begin(), inv_mass_trip.end());

    // Reset correction matrices containers.
    this->dirichlet_corr_mats_.clear(); this->dirichlet_corr_mats_.reserve(dirichlet_conds.size());
    this->loading_corr_mats_.clear(); this->loading_corr_mats_.reserve(loading_conds.size());

    // Set simplified to true until full EBCIEM is implemented.
    if (!is_simplified) {
        std::cout << Logger::Warning("Full EBCIEM not implemented yet. Simplified version will be used.") << std::endl;
        is_simplified = true;
    }

    // Get the correction matrices.
    if (is_simplified) { // Use Simplified EBCIEM (SEBCIEM)

        // Dirichlet nodes correction matrices.
        for (const auto &dirichlet : dirichlet_conds) {
            // Add dirichlet boundary correction matrix.
            this->dirichlet_corr_mats_.emplace_back(this->SebciemCorrMat(model.TetrahedralMesh().NodesNum(),
                                                                     dirichlet.NodesIds(), inv_mass));
        }

        // Loading nodes correction matrices.
        for (const auto &loading : loading_conds) {
            // Add loading boundary correction matrix.
            this->loading_corr_mats_.emplace_back(this->SebciemCorrMat(model.TetrahedralMesh().NodesNum(),
                                                                     loading.NodesIds(), inv_mass));
        }

    }
    else { // Use Full EBCIEM (SEBCIEM)
        //nothing to do for now.
    } // Get the correction matrices.


}


Eigen::SparseMatrix<double> Ebciem::SebciemCorrMat(const int &model_nodes_num, const std::vector<int> &bound_nodes_ids,
                                                   const Eigen::SparseMatrix<double> &model_inv_mass_mat) const
{
    // Initialize V matrix for EBCIEM.
    Eigen::SparseMatrix<double> v_mat(model_nodes_num, bound_nodes_ids.size());
    v_mat.setZero();

    // Iterate over boundary nodes.
    for (const auto &node_id : bound_nodes_ids) {
        auto it_id = &node_id - &bound_nodes_ids[0];
        v_mat.col(it_id) = this->nodal_mmls_.ShapeFunction().col(node_id);
    }

    // Compute decomposition of the v_mat_transpose*inv_mass*v_mat product.
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > decomp_prod;
    decomp_prod.compute(v_mat.transpose()*model_inv_mass_mat*v_mat);

    // Compute identity matrix.
    Eigen::SparseMatrix<double> identity(bound_nodes_ids.size(), bound_nodes_ids.size());
    identity.setIdentity();

    // Compute the inverse matrix of the product.
    Eigen::SparseMatrix<double> inv_prod = decomp_prod.solve(identity);

    // Compute the correction matrix.
    Eigen::SparseMatrix<double> corr_mat = model_inv_mass_mat*v_mat*inv_prod;

    // Return the correction matrix.
    return corr_mat;

}




} //end of namespace ExplicitSim
