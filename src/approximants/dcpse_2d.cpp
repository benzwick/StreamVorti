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

#include "StreamVorti/approximants/dcpse_2d.hpp"

#include <filesystem>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace StreamVorti {

void Dcpse2d::Update()
{
    // If something goes wrong we abort; but for now everything is OK so
    bool abort = false;

    mfem::StopWatch timer;
    timer.Start();
    std::cout << "DCPSE: update derivative matrices" << std::endl;

    mfem::GridFunction geom_nodes = this->SupportNodes();
    std::vector<std::vector<int> > support_nodes_ids = this->NeighborIndices();
    std::vector<double> support_radiuses = this->SupportRadiuses();

    const mfem::FiniteElementSpace *fes = geom_nodes.FESpace();
    int nnodes = fes->GetNDofs();

    // Initialize shape function derivative matrices.
    this->sh_func_dx_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dy_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dxx_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dyy_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dxy_ = mfem::SparseMatrix(nnodes, nnodes);

    // Keep track of min/max number of support nodes
    // TODO: Save these to paraview:
    // Add as members to DCPSE class as grid functions!
    // then let the user decide if they want to output them or not!
    // - number of neighbors
    // - cond(A)
    int min_supp_nodes = INT_MAX;
    int max_supp_nodes = INT_MIN;
    double min_cond_A1 = DBL_MAX;
    double max_cond_A1 = DBL_MIN;

    // Iterate over all the nodes of the grid.
    for (int node_id = 0; node_id < nnodes; ++node_id)
    {
        // for(const auto &node: geom_nodes) {
        //     auto node_id = &node - &geom_nodes[0];

        int nsupp = support_nodes_ids[node_id].size();
        min_supp_nodes = std::min(min_supp_nodes, nsupp);
        max_supp_nodes = std::max(max_supp_nodes, nsupp);
        if (nsupp < 1)
        {
            // std::cout << "Node " << node_id << " has no neighbors!" << std::endl;
            MFEM_WARNING("Node " << node_id << " has no neighbors. Abort!");
            abort = true;
            goto finish;
        }

        double node_X = geom_nodes(fes->DofToVDof(node_id, 0));
        double node_Y = geom_nodes(fes->DofToVDof(node_id, 1));

        // Monomial basis - Vandermonde matrices.
        Eigen::MatrixXd V1(nsupp, 6);
        Eigen::MatrixXd V2(nsupp, 5);

        Eigen::VectorXd expWd(nsupp);
        std::vector<Eigen::Triplet<double> > expW;
        expW.reserve(nsupp);

        double x = 0.;
        double y = 0.;
        double Wd = 0.;
        double val = 0.; double temp = 0.;
        double epsilon = 0.;

        // Iterate over node's neighbors.
        for (long unsigned int it = 0; it < nsupp; ++it)
        {
        // for (const auto &neigh_id : support_nodes_ids[node_id]) {
        //     auto it = &neigh_id - &support_nodes_ids[node_id][0];

            int neigh_id = support_nodes_ids[node_id][it];
            double neigh_X = geom_nodes(fes->DofToVDof(neigh_id, 0));
            double neigh_Y = geom_nodes(fes->DofToVDof(neigh_id, 1));

            // FIXME: should scale by 0.3 or not?
            epsilon = 0.3 * support_radiuses[node_id];
            x = (node_X - neigh_X) / epsilon;
            y = (node_Y - neigh_Y) / epsilon;

            //Set V1 vandermonde matrix for 1st derivatives.
            V1.coeffRef(it, 0) = 1.;
            V1.coeffRef(it, 1) = x;
            V1.coeffRef(it, 2) = y;
            V1.coeffRef(it, 3) = x*x;
            V1.coeffRef(it, 4) = x*y;
            V1.coeffRef(it, 5) = y*y;

            double neigh_dist_sq = (neigh_X - node_X) *
                                   (neigh_X - node_X) +
                                   (neigh_Y - node_Y) *
                                   (neigh_Y - node_Y);
            //Set expWd vector.
            Wd = - (neigh_dist_sq / (epsilon*epsilon) );
            expWd(it) = std::exp(Wd);

            //Set expW triplet.
            val = - (neigh_dist_sq / (2.*epsilon*epsilon) );
            temp = std::exp(val);
            expW.emplace_back(Eigen::Triplet<double>(it, it, temp) );

            //Set V2 vandermonde matrix for 2nd derivatives.
            V2(it, 0) = x;
            V2(it, 1) = y;
            V2(it, 2) = x*x;
            V2(it, 3) = x*y;
            V2(it, 4) = y*y;

        } // End Iterate over node's neighbors.

        Eigen::SparseMatrix<double> E(nsupp, nsupp);
        E.setFromTriplets(expW.begin(), expW.end());

        //1st order.
        Eigen::Matrix<double, 6, 6> A1;
        Eigen::MatrixXd B1(nsupp, 6);
        Eigen::MatrixXd B1_trans(6, nsupp);
        B1 = E*V1;
        B1_trans = B1.transpose();
        A1 = B1_trans*B1;

        // Condition number of 1st order A matrix
        {
            // double condA1 = Eigen::pseudoInverse(A1).norm() * A1.norm(); // not recommended?
            //
            // https://forum.kde.org/viewtopic.php?f=74&t=117430#p292018
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(A1);
            double condA1 = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
            min_cond_A1 = std::min(min_cond_A1, condA1);
            max_cond_A1 = std::max(max_cond_A1, condA1);
            // Check if cond(A) is within allowable limits
            if (condA1 > this->cond_A_limit_warn || condA1 <= 0.0)
            {
                // TODO: create string first then call the macro
                // if (num_condA_warnings > 3)
                //// NOTE: this is too noisy so comment out for now
                // MFEM_WARNING("cond(A1) = " << condA1
                //              << " at node " << node_id << ".");
                ////
                if (condA1 > this->cond_A_limit_abort || condA1 <= 0.0)
                {
                    MFEM_WARNING("cond(A1) exceeds limit of " << cond_A_limit_abort
                                 << " at node " << node_id << "."
                                 << " Abort!");
                    abort = true;
                    goto finish;
                }
            }
        }

        //x & y derivatives vectors.
        Eigen::Matrix<double, 6, 1> bx, by;
        bx.setZero(); bx(1) = -1.;
        by.setZero(); by(2) = -1.;

        //Solve linear systems.
        Eigen::Matrix<double, 6, 1> aTx, aTy;
        Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 6, 6> > sol1(A1);

        aTx = sol1.solve(bx);
        aTy = sol1.solve(by);

        Eigen::VectorXd coeffsx(nsupp);
        coeffsx = V1*aTx;
        coeffsx = coeffsx.cwiseProduct(expWd);

        Eigen::VectorXd coeffsy(nsupp);
        coeffsy = V1*aTy;
        coeffsy = coeffsy.cwiseProduct(expWd);

        //2nd order.
        Eigen::Matrix<double, 5, 5> A2;
        Eigen::MatrixXd B2(nsupp, 5);
        Eigen::MatrixXd B2_trans(5, nsupp);
        B2 = E*V2;
        B2_trans = B2.transpose();
        A2 = B2_trans*B2;

        //xx & yy derivatives vectors
        Eigen::Matrix<double, 5, 1> bxx, byy, bxy;
        bxx.setZero(); bxx(2) = 2.;
        byy.setZero(); byy(4) = 2.;
        bxy.setZero(); bxy(3) = 1.;

        //solve linear systems
        Eigen::VectorXd aTxx(5), aTyy(5), aTxy(5);
        Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 5, 5> > sol2(A2);

        aTxx = sol2.solve(bxx);
        aTyy = sol2.solve(byy);
        aTxy = sol2.solve(bxy);

        Eigen::VectorXd coeffsxx(nsupp);
        coeffsxx = V2*aTxx;
        coeffsxx = coeffsxx.cwiseProduct(expWd);

        Eigen::VectorXd coeffsyy(nsupp);
        coeffsyy = V2*aTyy;
        coeffsyy = coeffsyy.cwiseProduct(expWd);

        Eigen::VectorXd coeffsxy(nsupp);
        coeffsxy = V2*aTxy;
        coeffsxy = coeffsxy.cwiseProduct(expWd);

        //Shape functions derivatives components
        Eigen::VectorXd valX(nsupp);
        Eigen::VectorXd valY(nsupp);
        Eigen::VectorXd valXX(nsupp);
        Eigen::VectorXd valYY(nsupp);
        Eigen::VectorXd valXY(nsupp);

        valX.setZero(); valY.setZero();
        valXX.setZero(); valYY.setZero();
        valXY.setZero();

        double sumCoeffsx = 0., sumCoeffsy = 0.,
          sumCoeffsxx = 0., sumCoeffsyy = 0.,
          sumCoeffsxy = 0.;
        for(const auto &neigh_id : support_nodes_ids[node_id]) {
            auto it = &neigh_id - &support_nodes_ids[node_id][0];

            valX(it) = coeffsx(it) / epsilon;
            sumCoeffsx += coeffsx(it);

            valY(it) = coeffsy(it) / epsilon;
            sumCoeffsy += coeffsy(it);

            valXX(it) = coeffsxx(it) / (epsilon*epsilon);
            sumCoeffsxx += coeffsxx(it);

            valYY(it) = coeffsyy(it) / (epsilon*epsilon);
            sumCoeffsyy += coeffsyy(it);

            valXY(it) = coeffsxy(it) / (epsilon*epsilon);
            sumCoeffsxy += coeffsxy(it);
        }

        for(const auto &neigh_id : support_nodes_ids[node_id]) {
            auto it = &neigh_id - &support_nodes_ids[node_id][0];
            if(neigh_id == node_id) {
                valX(it) = valX(it) + (sumCoeffsx/epsilon);
                valY(it) = valY(it) + (sumCoeffsy/epsilon);
                valXX(it) = valXX(it) - (sumCoeffsxx/(epsilon*epsilon));
                valYY(it) = valYY(it) - (sumCoeffsyy/(epsilon*epsilon));
                valXY(it) = valXY(it) - (sumCoeffsxy/(epsilon*epsilon));
            }
        }

        //Store in triplets (center_node id, neighbour id, shFunc value).
        for(const auto &neigh_id : support_nodes_ids[node_id]) {
            auto it = &neigh_id - &support_nodes_ids[node_id][0];
            this->sh_func_dx_.Add(node_id, neigh_id, valX(it));
            this->sh_func_dy_.Add(node_id, neigh_id, valY(it));
            this->sh_func_dxx_.Add(node_id, neigh_id, valXX(it));
            this->sh_func_dyy_.Add(node_id, neigh_id, valYY(it));
            this->sh_func_dxy_.Add(node_id, neigh_id, valXY(it));
        }
    } // End Iterate over all the nodes of the grid.


    // Finalize shape function derivative matrices.
    this->sh_func_dx_.Finalize();
    this->sh_func_dy_.Finalize();
    this->sh_func_dxx_.Finalize();
    this->sh_func_dyy_.Finalize();
    this->sh_func_dxy_.Finalize();

finish:
    std::cout << "DCPSE: Min number of support nodes: " << min_supp_nodes << std::endl;
    std::cout << "DCPSE: Max number of support nodes: " << max_supp_nodes << std::endl;
    std::cout << "DCPSE: Min condition number of A1 matrix: " << min_cond_A1 << std::endl;
    std::cout << "DCPSE: Max condition number of A1 matrix: " << max_cond_A1 << std::endl;

    std::cout << "DCPSE: Execution time for DC PSE derivatives: "
              << timer.RealTime() << " s" << std::endl;

    if (abort)
    {
        MFEM_ABORT("DCPSE: Something bad happened.");
    }
}


void Dcpse2d::SaveDerivToFile(const std::string &deriv, const std::string &filename) const
{
    mfem::SparseMatrix derivative;

    if      (deriv == "dx")  { derivative = this->sh_func_dx_; }
    else if (deriv == "dy")  { derivative = this->sh_func_dy_; }
    else if (deriv == "dxx") { derivative = this->sh_func_dxx_; }
    else if (deriv == "dyy") { derivative = this->sh_func_dyy_; }
    else if (deriv == "dxy") { derivative = this->sh_func_dxy_; }

    if (derivative.Height() == 0)
    {
        MFEM_ABORT( "Could not save DCPSE derivative. "
                    "Derivatives have not been computed." );
    }

    //Initialize the path of the exporting file.
    std::string path = "";

    // Position of the last slash in the exporting file's name.
    std::size_t last_slash = filename.find_last_of("/\\");

    // Get the path directory of the exporting file name.
    if (last_slash != std::string::npos) {
        path = filename.substr(0, last_slash);
    }

    // Create the path's directory if it doesn't exist.
    std::filesystem::path dir(path);
    if (!path.empty() && !std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
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

    derivative.PrintMatlab(out);

    out.close();

}

} // namespace StreamVorti
