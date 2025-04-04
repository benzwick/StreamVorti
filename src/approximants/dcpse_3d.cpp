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

#include "StreamVorti/approximants/dcpse_3d.hpp"

#include <filesystem>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace StreamVorti {

void Dcpse3d::Update()
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
    this->sh_func_dz_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dxx_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dyy_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dzz_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dxy_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dxz_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dyz_ = mfem::SparseMatrix(nnodes, nnodes);

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
        double node_Z = geom_nodes(fes->DofToVDof(node_id, 2));

        // Monomial basis - Vandermonde matrices.
        Eigen::MatrixXd V1(nsupp, 10);
        Eigen::MatrixXd V2(nsupp, 9);

        Eigen::VectorXd expWd(nsupp);
        std::vector<Eigen::Triplet<double> > expW;
        expW.reserve(nsupp);

        double x = 0.;
        double y = 0.;
        double z = 0.;
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
            double neigh_Z = geom_nodes(fes->DofToVDof(neigh_id, 2));

            // FIXME: should scale by 0.3 or not?
            epsilon = 0.3 * support_radiuses[node_id];
            x = (node_X - neigh_X) / epsilon;
            y = (node_Y - neigh_Y) / epsilon;
            z = (node_Z - neigh_Z) / epsilon;

            //Set V1 vandermonde matrix for 1st derivatives.
            V1.coeffRef(it, 0) = 1.;
            V1.coeffRef(it, 1) = x;
            V1.coeffRef(it, 2) = y;
            V1.coeffRef(it, 3) = z;
            V1.coeffRef(it, 4) = x*x;
            V1.coeffRef(it, 5) = x*y;
            V1.coeffRef(it, 6) = x*z;
            V1.coeffRef(it, 7) = y*y;
            V1.coeffRef(it, 8) = y*z;
            V1.coeffRef(it, 9) = z*z;

            double neigh_dist_sq = (neigh_X - node_X) *
                                   (neigh_X - node_X) +
                                   (neigh_Y - node_Y) *
                                   (neigh_Y - node_Y) +
                                   (neigh_Z - node_Z) *
                                   (neigh_Z - node_Z);
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
            V2(it, 2) = z;
            V2(it, 3) = x*x;
            V2(it, 4) = x*y;
            V2(it, 5) = x*z;
            V2(it, 6) = y*y;
            V2(it, 7) = y*z;
            V2(it, 8) = z*z;

        } // End Iterate over node's neighbors.

        Eigen::SparseMatrix<double> E(nsupp, nsupp);
        E.setFromTriplets(expW.begin(), expW.end());

        //1st order.
        Eigen::Matrix<double, 10, 10> A1;
        Eigen::MatrixXd B1(nsupp, 10);
        Eigen::MatrixXd B1_trans(10, nsupp);
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
                MFEM_WARNING("cond(A1) = " << condA1
                             << " at node " << node_id << ".");
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
        Eigen::Matrix<double, 10, 1> bx, by, bz;
        bx.setZero(); bx(1) = -1.;
        by.setZero(); by(2) = -1.;
        bz.setZero(); bz(3) = -1.;

        //Solve linear systems.
        Eigen::Matrix<double, 10, 1> aTx, aTy, aTz;
        Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 10, 10> > sol1(A1);

        aTx = sol1.solve(bx);
        aTy = sol1.solve(by);
        aTz = sol1.solve(bz);

        Eigen::VectorXd coeffsx(nsupp);
        coeffsx = V1*aTx;
        coeffsx = coeffsx.cwiseProduct(expWd);

        Eigen::VectorXd coeffsy(nsupp);
        coeffsy = V1*aTy;
        coeffsy = coeffsy.cwiseProduct(expWd);

        Eigen::VectorXd coeffsz(nsupp);
        coeffsz = V1*aTz;
        coeffsz = coeffsz.cwiseProduct(expWd);

        //2nd order.
        Eigen::Matrix<double, 9, 9> A2;
        Eigen::MatrixXd B2(nsupp, 9);
        Eigen::MatrixXd B2_trans(9, nsupp);
        B2 = E*V2;
        B2_trans = B2.transpose();
        A2 = B2_trans*B2;

        //xx & yy derivatives vectors
        Eigen::Matrix<double, 9, 1> bxx, byy, bzz, bxy, bxz, byz;
        bxx.setZero(); bxx(3) = 2.;
        byy.setZero(); byy(6) = 2.;
        bzz.setZero(); bzz(8) = 2.;
        bxy.setZero(); bxy(4) = 1.;
        bxz.setZero(); bxz(5) = 1.;
        byz.setZero(); byz(7) = 1.;

        //solve linear systems
        Eigen::VectorXd aTxx(9), aTyy(9), aTzz(9), aTxy(9), aTxz(9), aTyz(9);
        Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 9, 9> > sol2(A2);

        aTxx = sol2.solve(bxx);
        aTyy = sol2.solve(byy);
        aTzz = sol2.solve(bzz);
        aTxy = sol2.solve(bxy);
        aTxz = sol2.solve(bxz);
        aTyz = sol2.solve(byz);

        Eigen::VectorXd coeffsxx(nsupp);
        coeffsxx = V2*aTxx;
        coeffsxx = coeffsxx.cwiseProduct(expWd);

        Eigen::VectorXd coeffsyy(nsupp);
        coeffsyy = V2*aTyy;
        coeffsyy = coeffsyy.cwiseProduct(expWd);

        Eigen::VectorXd coeffszz(nsupp);
        coeffszz = V2*aTzz;
        coeffszz = coeffszz.cwiseProduct(expWd);

        Eigen::VectorXd coeffsxy(nsupp);
        coeffsxy = V2*aTxy;
        coeffsxy = coeffsxy.cwiseProduct(expWd);

        Eigen::VectorXd coeffsxz(nsupp);
        coeffsxz = V2*aTxz;
        coeffsxz = coeffsxz.cwiseProduct(expWd);

        Eigen::VectorXd coeffsyz(nsupp);
        coeffsyz = V2*aTyz;
        coeffsyz = coeffsyz.cwiseProduct(expWd);

        //Shape functions derivatives components
        Eigen::VectorXd valX(nsupp);
        Eigen::VectorXd valY(nsupp);
        Eigen::VectorXd valZ(nsupp);
        Eigen::VectorXd valXX(nsupp);
        Eigen::VectorXd valYY(nsupp);
        Eigen::VectorXd valZZ(nsupp);
        Eigen::VectorXd valXY(nsupp);
        Eigen::VectorXd valXZ(nsupp);
        Eigen::VectorXd valYZ(nsupp);

        valX.setZero(); valY.setZero(); valZ.setZero();
        valXX.setZero(); valYY.setZero(); valZZ.setZero();
        valXY.setZero(); valXZ.setZero(); valYZ.setZero();

        double sumCoeffsx = 0., sumCoeffsy = 0., sumCoeffsz = 0.,
          sumCoeffsxx = 0., sumCoeffsyy = 0., sumCoeffszz = 0.,
          sumCoeffsxy = 0., sumCoeffsxz = 0., sumCoeffsyz = 0.;
        for(const auto &neigh_id : support_nodes_ids[node_id]) {
            auto it = &neigh_id - &support_nodes_ids[node_id][0];

            valX(it) = coeffsx(it) / epsilon;
            sumCoeffsx += coeffsx(it);

            valY(it) = coeffsy(it) / epsilon;
            sumCoeffsy += coeffsy(it);

            valZ(it) = coeffsz(it) / epsilon;
            sumCoeffsz += coeffsz(it);

            valXX(it) = coeffsxx(it) / (epsilon*epsilon);
            sumCoeffsxx += coeffsxx(it);

            valYY(it) = coeffsyy(it) / (epsilon*epsilon);
            sumCoeffsyy += coeffsyy(it);

            valZZ(it) = coeffszz(it) / (epsilon*epsilon);
            sumCoeffszz += coeffszz(it);

            valXY(it) = coeffsxy(it) / (epsilon*epsilon);
            sumCoeffsxy += coeffsxy(it);

            valXZ(it) = coeffsxz(it) / (epsilon*epsilon);
            sumCoeffsxz += coeffsxz(it);

            valYZ(it) = coeffsyz(it) / (epsilon*epsilon);
            sumCoeffsyz += coeffsyz(it);
        }

        for(const auto &neigh_id : support_nodes_ids[node_id]) {
            auto it = &neigh_id - &support_nodes_ids[node_id][0];
            if(neigh_id == node_id) {
                valX(it) = valX(it) + (sumCoeffsx/epsilon);
                valY(it) = valY(it) + (sumCoeffsy/epsilon);
                valZ(it) = valZ(it) + (sumCoeffsy/epsilon);
                valXX(it) = valXX(it) - (sumCoeffsxx/(epsilon*epsilon));
                valYY(it) = valYY(it) - (sumCoeffsyy/(epsilon*epsilon));
                valZZ(it) = valZZ(it) - (sumCoeffszz/(epsilon*epsilon));
                valXY(it) = valXY(it) - (sumCoeffsxy/(epsilon*epsilon));
                valXZ(it) = valXZ(it) - (sumCoeffsxz/(epsilon*epsilon));
                valYZ(it) = valYZ(it) - (sumCoeffsyz/(epsilon*epsilon));
            }
        }

        //Store in triplets (center_node id, neighbour id, shFunc value).
        for(const auto &neigh_id : support_nodes_ids[node_id]) {
            auto it = &neigh_id - &support_nodes_ids[node_id][0];
            this->sh_func_dx_.Add(node_id, neigh_id, valX(it));
            this->sh_func_dy_.Add(node_id, neigh_id, valY(it));
            this->sh_func_dz_.Add(node_id, neigh_id, valZ(it));
            this->sh_func_dxx_.Add(node_id, neigh_id, valXX(it));
            this->sh_func_dyy_.Add(node_id, neigh_id, valYY(it));
            this->sh_func_dzz_.Add(node_id, neigh_id, valZZ(it));
            this->sh_func_dxy_.Add(node_id, neigh_id, valXY(it));
            this->sh_func_dxz_.Add(node_id, neigh_id, valXZ(it));
            this->sh_func_dyz_.Add(node_id, neigh_id, valYZ(it));
        }
    } // End Iterate over all the nodes of the grid.


    // Finalize shape function derivative matrices.
    this->sh_func_dx_.Finalize();
    this->sh_func_dy_.Finalize();
    this->sh_func_dz_.Finalize();
    this->sh_func_dxx_.Finalize();
    this->sh_func_dyy_.Finalize();
    this->sh_func_dzz_.Finalize();
    this->sh_func_dxy_.Finalize();
    this->sh_func_dxz_.Finalize();
    this->sh_func_dyz_.Finalize();

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


void Dcpse3d::SaveDerivToFile(const std::string &deriv, const std::string &filename) const
{
    mfem::SparseMatrix derivative;

    if      (deriv == "dx")  { derivative = this->sh_func_dx_; }
    else if (deriv == "dy")  { derivative = this->sh_func_dy_; }
    else if (deriv == "dz")  { derivative = this->sh_func_dz_; }
    else if (deriv == "dxx") { derivative = this->sh_func_dxx_; }
    else if (deriv == "dyy") { derivative = this->sh_func_dyy_; }
    else if (deriv == "dzz") { derivative = this->sh_func_dzz_; }
    else if (deriv == "dxy") { derivative = this->sh_func_dxy_; }
    else if (deriv == "dxz") { derivative = this->sh_func_dxz_; }
    else if (deriv == "dyz") { derivative = this->sh_func_dyz_; }

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
