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


#include "StreamVorti/approximants/dcpse_2d.hpp"


namespace StreamVorti {


Dcpse2d::Dcpse2d()
{}


Dcpse2d::~Dcpse2d()
{}


void Dcpse2d::ComputeDerivs(const std::vector<Node> &geom_nodes,
                            const std::vector<std::vector<int> > &support_nodes_ids,
                            const std::vector<double> &support_radiuses)
{

    // Initialize shape function derivative matrices.
    this->sh_func_dx_ = mfem::SparseMatrix(geom_nodes.size(), geom_nodes.size());
    this->sh_func_dy_ = mfem::SparseMatrix(geom_nodes.size(), geom_nodes.size());
    this->sh_func_dxx_ = mfem::SparseMatrix(geom_nodes.size(), geom_nodes.size());
    this->sh_func_dyy_ = mfem::SparseMatrix(geom_nodes.size(), geom_nodes.size());
    this->sh_func_dxy_ = mfem::SparseMatrix(geom_nodes.size(), geom_nodes.size());

    // Iterate over all the nodes of the grid.
    for(const auto &node: geom_nodes) {
        auto node_id = &node - &geom_nodes[0];

        // Monomial basis - Vandermonde matrices.
        Eigen::MatrixXd V1(support_nodes_ids[node_id].size(),6);
        Eigen::MatrixXd V2(support_nodes_ids[node_id].size(),5);

        Eigen::VectorXd expWd(support_nodes_ids[node_id].size());
        std::vector<Eigen::Triplet<double> > expW;
        expW.reserve(support_nodes_ids[node_id].size());

        double x = 0.; double y = 0.;
        double Wd = 0.;
        double val = 0.; double temp = 0.;
        double epsilon = 0.;

        // Iterate over node's neighbors.
        for (const auto &neigh_id : support_nodes_ids[node_id]) {
            auto it = &neigh_id - &support_nodes_ids[node_id][0];

            epsilon = support_radiuses[node_id];
            x = (node.Coordinates().X() - geom_nodes[neigh_id].Coordinates().X()) / epsilon;
            y = (node.Coordinates().Y() - geom_nodes[neigh_id].Coordinates().Y()) / epsilon;

            //Set V1 vandermonde matrix for 1st derivatives.
            V1.coeffRef(it, 0) = 1.;
            V1.coeffRef(it, 1) = x;
            V1.coeffRef(it, 2) = y;
            V1.coeffRef(it, 3) = x*x;
            V1.coeffRef(it, 4) = x*y;
            V1.coeffRef(it, 5) = y*y;

            double neigh_dist_sq = (geom_nodes[neigh_id].Coordinates().X() - node.Coordinates().X()) *
                                   (geom_nodes[neigh_id].Coordinates().X() - node.Coordinates().X()) +
                                   (geom_nodes[neigh_id].Coordinates().Y() - node.Coordinates().Y()) *
                                   (geom_nodes[neigh_id].Coordinates().Y() - node.Coordinates().Y());
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

        Eigen::SparseMatrix<double> E(support_nodes_ids[node_id].size(), support_nodes_ids[node_id].size());
        E.setFromTriplets(expW.begin(), expW.end());

        //1st order.
        Eigen::Matrix<double, 6, 6> A1;
        Eigen::MatrixXd B1(support_nodes_ids[node_id].size(), 6);
        Eigen::MatrixXd B1_trans(6, support_nodes_ids[node_id].size());
        B1 = E*V1;
        B1_trans = B1.transpose();
        A1 = B1_trans*B1;

        //x & y derivatives vectors.
        Eigen::Matrix<double, 6, 1> bx, by;
        bx(0) = 0.; bx(1) = -1.; bx(2) = 0.; bx(3) = 0.; bx(4)=0., bx(5)=0.;
        by(0) = 0.; by(1) = 0.; by(2) = -1.; by(3) = 0.; by(4)=0., by(5)=0.;

        //Solve linear systems.
        Eigen::Matrix<double, 6, 1> aTx, aTy;
        Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 6, 6> > sol1(A1);

        aTx = sol1.solve(bx);
        aTy = sol1.solve(by);

        Eigen::VectorXd coeffsx(support_nodes_ids[node_id].size());
        coeffsx = V1*aTx;
        coeffsx = coeffsx.cwiseProduct(expWd);

        Eigen::VectorXd coeffsy(support_nodes_ids[node_id].size());
        coeffsy = V1*aTy;
        coeffsy = coeffsy.cwiseProduct(expWd);

        //2nd order.
        Eigen::Matrix<double, 5, 5> A2;
        Eigen::MatrixXd B2(support_nodes_ids[node_id].size(), 5);
        Eigen::MatrixXd B2_trans(5, support_nodes_ids[node_id].size());
        B2 = E*V2;
        B2_trans = B2.transpose();
        A2 = B2_trans*B2;

        //xx & yy derivatives vectors
        Eigen::Matrix<double, 5, 1> bxx, byy, bxy;
        bxx(0) = 0.; bxx(1) = 0.; bxx(2) = 2.; bxx(3) = 0.; bxx(4)=0.;
        byy(0) = 0.; byy(1) = 0.; byy(2) = 0.; byy(3) = 0.; byy(4)=2.;
        bxy(0) = 0.; bxy(1) = 0.; bxy(2) = 0.; bxy(3) = 1.; bxy(4)=0.;

        //solve linear systems
        Eigen::VectorXd aTxx(5), aTyy(5), aTxy(5);
        Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 5, 5> > sol2(A2);

        aTxx = sol2.solve(bxx);
        aTyy = sol2.solve(byy);
        aTxy = sol2.solve(bxy);

        Eigen::VectorXd coeffsxx(support_nodes_ids[node_id].size());
        coeffsxx = V2*aTxx;
        coeffsxx = coeffsxx.cwiseProduct(expWd);

        Eigen::VectorXd coeffsyy(support_nodes_ids[node_id].size());
        coeffsyy = V2*aTyy;
        coeffsyy = coeffsyy.cwiseProduct(expWd);

        Eigen::VectorXd coeffsxy(support_nodes_ids[node_id].size());
        coeffsxy = V2*aTxy;
        coeffsxy = coeffsxy.cwiseProduct(expWd);

        //Shape functions derivatives components
        Eigen::VectorXd valX(support_nodes_ids[node_id].size());
        Eigen::VectorXd valY(support_nodes_ids[node_id].size());
        Eigen::VectorXd valXX(support_nodes_ids[node_id].size());
        Eigen::VectorXd valYY(support_nodes_ids[node_id].size());
        Eigen::VectorXd valXY(support_nodes_ids[node_id].size());

        valX.setZero(); valY.setZero(); valXX.setZero();
        valYY.setZero(); valXY.setZero();

        double sumCoeffsx = 0., sumCoeffsy = 0., sumCoeffsxx = 0., sumCoeffsyy = 0., sumCoeffsxy = 0.;
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

}


void Dcpse2d::SaveDerivToFile(const std::string &deriv, const std::string &filename) const
{
    mfem::SparseMatrix derivative;

    if (deriv == "dx") { derivative = this->sh_func_dx_; }
    else if (deriv == "dy") { derivative = this->sh_func_dy_; }
    else if (deriv == "dxx") { derivative = this->sh_func_dxx_; }
    else if (deriv == "dyy") { derivative = this->sh_func_dyy_; }
    else if (deriv == "dxy") { derivative = this->sh_func_dxy_; }


    if (derivative.Height() == 0) {
        throw std::runtime_error(Logger::Error("Could not save DCPSE derivative. "
                                               "Derivatives have not been computed.").c_str());
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

    derivative.PrintMatlab(out);

    out.close();

}


} //end of namespace StreamVorti
