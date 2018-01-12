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
 *      Konstantinos A. MOUNTRIS
 */

#include <StreamVorti/stream_vorti.hpp>

#include <boost/filesystem.hpp>

#include <cstddef>
#include <string>
#include <fstream>

using namespace StreamVorti;

int main(int argc, char *argv[]) {

//    try {

//        // Initialize configuration manager.
//        ConfigManager *config = new ConfigManager();

//        // Initialize configuration filename empty to be read during the ExplicitSim execution.
//        std::string config_filename = "";

//        // Check if configuration file was provided during execution.
//        if (argc == 1) {
//            std::cout << Logger::Warning("No input file was specified. Type configuration filename with absolute path.\n"
//                         "Otherwise tap '-g' to generate sample configuration file or '-q' to exit ExplicitSim.\nInput filename: ");

//            // Read configuration file given by the user.
//            std::cin >> config_filename;

//            if (config_filename == "-g") {
//                std::cout << Logger::Message("Give text file name [.ini] with absolute path to store the sample configuration "
//                             "file\nor tap '-t' to print in terminal.\nSample filename: ");

//                std::string sample_filename = "";

//                std::cin >> sample_filename;

//                if (sample_filename == "-t") {
//                    std::cout << config->PrintSampleFile() << std::endl;
//                    std::cout << Logger::Message("Save the sample configuration in a text file [.ini], edit it according to your simulation,"
//                                 " and relaunch ExplicitSimRun passing your configuration file as argument.\n");
//                    return EXIT_SUCCESS;
//                }
//                else {
//                    // Initialize the path of the sample configuration file.
//                    std::string sample_path = "";

//                    // Position of the last slash in the sample configuration file.
//                    std::size_t last_slash = sample_filename.find_last_of("/\\");

//                    // Get the path directory of the sample configuration file.
//                    if (last_slash != std::string::npos) {
//                        sample_path = sample_filename.substr(0, last_slash);
//                    }

//                    // Create the path's directory if it doesn't exist.
//                    boost::filesystem::path dir(sample_path);
//                    if (!sample_path.empty() && !boost::filesystem::exists(dir)) {
//                        boost::filesystem::create_directories(dir);
//                    }

//                    // Position of the last slash in the exporting file's name.
//                    std::size_t last_dot = sample_filename.find_last_of(".");

//                    // Initialize sample configuration file extension.
//                    std::string sample_ext = "";

//                    // Get sample configuration file extension.
//                    if (last_dot != std::string::npos) {
//                        sample_ext = sample_filename.substr(last_dot);
//                    }

//                    // Add extension if necessary
//                    if (sample_ext != ".ini") { sample_filename += ".ini"; }

//                    // Output the sample file.
//                    std::ofstream sample_output(sample_filename, std::ios_base::out | std::ios_base::trunc);
//                    sample_output << config->PrintSampleFile();
//                    std::cout << Logger::Message("Sample configuration file saved at: ") << sample_filename << std::endl;
//                    std::cout << Logger::Message("Edit the sample configuration according to your simulation and "
//                                 "relaunch ExplicitSimRun passing your configuration file as argument.\n");
//                    return EXIT_SUCCESS;
//                }

//            }

//            if (config_filename == "-q") { std::cout << Logger::Message("User requested termination. See you soon!\n"); exit(0); }
//            std::cout << config_filename << std::endl;
//            exit(0);
//        }
//        else { config_filename = argv[1]; }

//        // Read configuration file.
//        config->ReadConfigFile(config_filename);

//        std::cout << "<< Welcome to ExplicitSim >>\n";
//        std::cout << Logger::Message("Loading configuration file: ") << config_filename << std::endl;

//        // Profiling spent time in ExplicitSim
//        Timer timer;

//        // Set the weak form model.
//        WeakModel3D model;
//        model.LoadMeshRepresentation(config->RetrieveArgument<std::string>("Model.MeshFile"));
//        model.CreateGridRepresentation();

//        std::cout << Logger::Message("Model has nodes: ") << model.TetrahedralMesh().Nodes().size() <<
//                     " and elements: " << model.TetrahedralMesh().Elements().size() << std::endl;

//        // Set support domain.
//        SupportDomain support;
//        support.SetInfluenceNodes(model.TetrahedralMesh().Nodes());
//        support.SetInfluenceTetrahedra(model.TetrahedralMesh().Elements());

//        // Compute the influnce radiuses of the support domain.
//        support.ComputeInfluenceNodesRadiuses(config->RetrieveArgument<double>("ShapeFunction.DilatationCoefficient"));

//        // Setting options for integration points generation.
//        IntegOptions options;
//        options.is_adaptive_ = config->RetrieveArgument<bool>("IntegrationOptions.Adaptive");
//        options.adaptive_eps_ = config->RetrieveArgument<double>("IntegrationOptions.AdaptiveEps");
//        options.adaptive_level_ = config->RetrieveArgument<int>("IntegrationOptions.AdaptiveLevel");
//        options.tetra_divisions_ = config->RetrieveArgument<int>("IntegrationOptions.TetrahedronDivisions");
//        options.integ_points_per_tetra_ = config->RetrieveArgument<int>("IntegrationOptions.IntegPointsPerTetrahedron");

//        if (options.is_adaptive_) {
//            std::cout << Logger::Message("Adaptive integration: ON\n");
//            std::cout << Logger::Message("Tetrahedron divisions used: ") << options.tetra_divisions_ << std::endl;
//            std::cout << Logger::Message("Number of integration points per tetrahedron division used: 4\n");
//        }
//        else {
//            std::cout << Logger::Message("Adaptive integration: OFF\n");
//            std::cout << Logger::Message("Number of integration points per element used: ")
//                      << options.integ_points_per_tetra_ << std::endl;
//        }

//        // Time only for integration points generation.
//        timer.Reset();

//        // Create integration points for the model.
//        model.CreateIntegrationPoints(options, support);
//        std::cout << Logger::Message("Model has integration points: ") << model.IntegrationPoints().PointsNum() << "\n";
//        std::cout << Logger::Message("Execution time for integration points generation: ") << timer.PrintElapsedTime() << "\n";

//        // Time only for closest nodes.
//        timer.Reset();
//        // Find influence nodes indices of integration points.
//        auto neighs_ids = support.CgalClosestNodesIdsTo(model.IntegrationPoints().Coordinates());

//        std::cout << Logger::Message("Execution time for neighbor nodes computation: ") << timer.PrintElapsedTime() << "\n";
//        std::cout << Logger::Message("The minimum and maximum number of support nodes: ")
//                  << support.MinSupportNodesIn(neighs_ids) << " - " << support.MaxSupportNodesIn(neighs_ids) << std::endl;

//        // Assign a neo-hookean material to the model.
//        NeoHookean material;
//        material.SetPointsNumber(model.IntegrationPoints().Coordinates().size());

//        // Set material parameters.
//        material.SetDensity(config->RetrieveArgument<double>("Material.Density"));
//        material.SetYoungModulus(config->RetrieveArgument<double>("Material.YoungModulus"));
//        material.SetPoissonRatio(config->RetrieveArgument<double>("Material.PoissonRatio"));

//        // Compute elastic parameters and wave speed.
//        material.ComputeLameLambdaMu();
//        material.ComputeBulkModulus();
//        material.ComputeWaveSpeed();

//        // Time only for shape functions.
//        timer.Reset();

//        // Shape Functions
//        Mmls3d mmls3d;
//        mmls3d.SetBasisFunctionType(config->RetrieveArgument<std::string>("ShapeFunction.BasisFunctionType"));
//        mmls3d.SetExactDerivativesMode(config->RetrieveArgument<bool>("ShapeFunction.UseExactDerivatives"));
//        mmls3d.ComputeShFuncAndDerivs(model.TetrahedralMesh().Nodes(), model.IntegrationPoints().Coordinates(),
//                                          neighs_ids, support.InfluenceNodesRadiuses());

//        std::cout << Logger::Message("Execution time for shape function generation: ") << timer.PrintElapsedTime() << "\n";

//        // The dynamic relaxation properties.
//        DynRelaxProp dr;
//        dr.SetEquilibriumTime(config->RetrieveArgument<double>("DynamicRelaxation.EquilibriumTime"));
//        dr.SetLoadConvRate(config->RetrieveArgument<double>("DynamicRelaxation.LoadConvRate"));
//        dr.SetAfterLoadConvRate(config->RetrieveArgument<double>("DynamicRelaxation.AfterLoadConvRate"));
//        dr.SetStopUpdateConvRateStepsNum(config->RetrieveArgument<int>("DynamicRelaxation.StopUpdateConvRateStepsNum"));
//        dr.SetConvRateDeviation(config->RetrieveArgument<double>("DynamicRelaxation.ConvRateDeviation"));
//        dr.SetForceDispUpdateStepsNum(config->RetrieveArgument<int>("DynamicRelaxation.ForceDispUpdateStepsNum"));
//        dr.SetStableConvRateStepsNum(config->RetrieveArgument<int>("DynamicRelaxation.StableConvRateStepsNum"));
//        dr.SetConvRateStopDeviation(config->RetrieveArgument<double>("DynamicRelaxation.ConvRateStopDeviation"));
//        dr.SetStopConvRateError(config->RetrieveArgument<double>("DynamicRelaxation.StopConvRateError"));
//        dr.SetStopAbsError(config->RetrieveArgument<double>("DynamicRelaxation.StopAbsError"));
//        dr.SetStopStepsNum(config->RetrieveArgument<int>("DynamicRelaxation.StopStepsNum"));

//        // Solve the model explicitly with MTLED.
//        Mtled solver;

//        // Compute the time steps.
//        solver.ComputeTimeSteps(material.WaveSpeed(), neighs_ids, mmls3d);

//        // Compute the mass of the model.
//        model.ComputeMass(material.Density(), solver.TimeSteps(), solver.MaxStep(),
//                          neighs_ids, config->RetrieveArgument<bool>("Model.MassScaling"));

//        // Compute the stable time step for the explicit solution.
//        if (config->RetrieveArgument<bool>("MTLED.UsePredefinedStableTimeStep")) {
//            solver.SetStableStep(config->RetrieveArgument<double>("MTLED.StableTimeStep"));
//        }
//        else { solver.ComputeStableStep(model.Mass(), model.IsMassScaled(), 1.5); }

//        // Set steps for progress save.
//        solver.SetSaveProgressSteps(config->RetrieveArgument<int>("MTLED.SaveProgressSteps"));

//        std::cout << Logger::Message("Minimum time step: ") << solver.MinStep() << " s\n";
//        std::cout << Logger::Message("Maximum time step: ") << solver.MaxStep() << " s\n";
//        std::cout << Logger::Message("Stable time step: ") << solver.StableStep() << " s\n";

//        // Set the loading condition.
//        LoadCurve load;
//        load.SetLoadTime(config->RetrieveArgument<double>("Loading.LoadTime"));
//        load.SetMaxDisplacement(config->RetrieveArgument<double>("Loading.MaxDisplacement"));
//        load.ComputeLoadStepsNum(solver.StableStep());
//        load.ComputeLoadStepDisplacements(solver.StableStep());
//        std::cout << Logger::Message("The loading steps number: ") << load.LoadStepsNum() << std::endl;

//        // Compute dynamic relaxation steps number for equilibrium.
//        dr.ComputeEquilibriumStepsNum(solver.StableStep());
//        std::cout << Logger::Message("The dynamic relaxation equilibrium steps number: ") << dr.EquilibriumStepsNum() << std::endl;

//        // Compute total time steps for solution.
//        solver.ComputeTotalTimeStepsNum(load.LoadStepsNum(), dr.EquilibriumStepsNum());
//        std::cout << Logger::Message("The total explicit solution time steps number: ") << solver.TotalTimeStepsNum() << std::endl;

//        // Time only for conditions imposition.
//        timer.Reset();

//        // Initialize the boundary and loading conditions handler.
//        ConditionsHandler cond_handler;

//        // Set fixed displacement conditions.
//        if (config->VarMap().count("Boundary.FixedName")) {
//            int fixed_name_num = config->OptionsNumInList<std::string>("Boundary.FixedName");
//            int fixed_axes_num = config->OptionsNumInList<std::string>("Boundary.FixedAxes");

//            // Check if same number of fixed displacement conditions names and fixed axes triplets are given.
//            if (fixed_name_num != fixed_axes_num) {
//                throw std::invalid_argument(Logger::Error("Number of given fixed displacements conditions does not "
//                                                          "match to given axes fixation triplets. "
//                                                          "Check the given configuration options.").c_str());
//            }

//            // Fixed axes conditionals.
//            bool is_x_fixed = false; bool is_y_fixed = false; bool is_z_fixed = false;

//            // Add fixed displacement conditions to the conditions handler.
//            for (int cond_num = 0; cond_num != fixed_axes_num; ++cond_num) {
//                // String stream to parse fixation conditions from string argument.
//                std::stringstream ss(config->RetrieveArgumentFromList<std::string>("Boundary.FixedAxes", cond_num));

//                if (!(ss >> is_x_fixed >> is_y_fixed >> is_z_fixed)) {
//                    std::string error = Logger::Error("Could not process fixed axes triplet "
//                                        "for condition:") + std::to_string(cond_num);
//                    throw std::invalid_argument(error.c_str());
//                }

//                cond_handler.AddDirichlet(is_x_fixed, is_y_fixed, is_z_fixed,
//                                          config->RetrieveArgumentFromList<std::string>("Boundary.FixedName", cond_num));

//            }
//        } // End of Set fixed displacement conditions.


//        // Add loading condition in conditions handler.
//        // Loading axes conditionals.
//        bool is_x_loaded = false; bool is_y_loaded = false; bool is_z_loaded = false;

//        // String stream to parse loading conditions from string argument.
//        std::stringstream ss(config->RetrieveArgument<std::string>("Loading.LoadAxes"));

//        if (!(ss >> is_x_loaded >> is_y_loaded >> is_z_loaded)) {
//            throw std::invalid_argument(Logger::Error("Could not process loading axes triplet.").c_str());
//        }

//        cond_handler.AddLoading(load, is_x_loaded, is_y_loaded, is_z_loaded,
//                                config->RetrieveArgument<std::string>("Loading.LoadName"));

//        // Extact nodes ids where conditions are imposed.
//        cond_handler.ExtractBoundaryNodeIds(model.TetrahedralMesh().NodeSets());

//        // Add EBCIEM correction in the conditions' handler.
//        if (config->RetrieveArgument<bool>("EBCIEM.UseEBCIEM")) {
//            std::cout << Logger::Message("EBCIEM Application: ON") << std::endl;
//            cond_handler.AddEbciem(model, support, mmls3d.BaseFunctionType(), mmls3d.UseExactDerivatives(),
//                                   config->RetrieveArgument<bool>("EBCIEM.UseSimplifiedVersion"));
//        }

//        std::cout << Logger::Message("Execution time for conditions initialization: ") << timer.PrintElapsedTime() << "\n";

//        // Solving with explicit dynamics.
//        std::cout << Logger::Message("Starting the explicit solution of the model...") << std::endl;
//        timer.Reset();

//        solver.Solve(model, neighs_ids, cond_handler, mmls3d, material, dr,
//                     config->RetrieveArgument<bool>("EBCIEM.UseEBCIEM"));
//        std::cout << Logger::Message("Execution time for explicit solution: ") << timer.PrintElapsedTime() << "\n";

//        // Apply shape function correction to the nodal displacements.
//        std::cout << Logger::Message("Shape function values application on nodal displacements started...") << std::endl;
//        timer.Reset();

//        // Find influence nodes inidices the nodes of the model's geometry.
//        auto nodal_neigh_ids = support.CgalClosestNodesIdsTo(model.TetrahedralMesh().NodeCoordinates());

//        // Compute the mmls approximants for the nodes of the model's geometry.
//        Mmls3d nodal_mmls;
//        nodal_mmls.SetBasisFunctionType(mmls3d.BaseFunctionType());
//        nodal_mmls.SetExactDerivativesMode(mmls3d.UseExactDerivatives());
//        nodal_mmls.ComputeShFuncAndDerivs(model.TetrahedralMesh().Nodes(),
//                                          model.TetrahedralMesh().NodeCoordinates(),
//                                          nodal_neigh_ids, support.InfluenceNodesRadiuses());

//        // Application of the shape functions on the displacements.
//        solver.ApplyShapeFuncToDisplacements(model, nodal_mmls);

//        std::cout << Logger::Message("Shape function values application on nodal displacements completed in: ")
//                  << timer.PrintElapsedTime() << "\n";

//        // Create the paraview output.
//        if (config->VarMap().count("Output.FilePath") && config->VarMap().count("Output.FileName")) {

//            timer.Reset();
//            std::cout << Logger::Message("Saving results...\n");

//            // Model output filename.
//            std::string output_filepath = config->RetrieveArgument<std::string>("Output.FilePath");
//            std::string output_filename = config->RetrieveArgument<std::string>("Output.FileName");

//            // Strip extension from output filename.
//            std::string out_ext = output_filename.substr(output_filename.find_last_of("."));
//            if (out_ext == ".vtu") { output_filename.erase(output_filename.end()-4, output_filename.end()); }

//            ParaviewExporter exporter;
//            exporter.CreateVtu(model);

//            // Export a .vtu file for each saved time step.
//            for (std::size_t i = 0; i != solver.SavedDisplacements().size(); ++i) {
//                exporter.AddVectorField(solver.SavedDisplacements()[i], "Displacements");
//                exporter.AddVectorField(solver.SavedForces()[i], "Forces");
//                exporter.Export(output_filepath+output_filename+std::to_string(i)+".vtu");
//                exporter.ClearVectorFields();
//                std::cout << Logger::Message("Saved ") << i+1 << "/" << solver.SavedDisplacements().size()
//                          << " stored model states at: " << output_filepath+output_filename+std::to_string(i)+".vtu\n";
//            }

//            if (config->VarMap().count("Output.AnimationName")) {
//                // Create animation.
//                std::string animation_filename = output_filepath + config->RetrieveArgument<std::string>("Output.AnimationName");
//                exporter.CreatePvdAnimation(output_filepath, solver.SavedDisplacements().size(), animation_filename);
//                std::cout << Logger::Message("Saved animation of the stored model states at: ") << animation_filename << std::endl;
//            }

//            std::cout << Logger::Message("Time for output: ") << timer.PrintElapsedTime() << "\n";
//        }

//        std::cout << Logger::Message("Simulation terminated successfully.") << std::endl;

//        // Release memory from configuration manager.
//        delete config;

//    }
//    catch (const std::invalid_argument &e) {
//        std::cerr << e.what() << std::endl;
//    }
//    catch (const std::runtime_error &e) {
//        std::cerr << e.what() << std::endl;
//    }
//    catch (const std::out_of_range &e) {
//       std::cerr << e.what() << std::endl;
//    }
//    catch (const std::bad_alloc &e) {
//       std::cerr << e.what() << std::endl;
//    }
//    catch (const boost::program_options::error &e) {
//        std::cerr << "[Boost program_options error] " << e.what() << std::endl;
//    }
//    catch (...) {
//        std::cerr << "[ExplicitSim Unknown exception]" << std::endl;
//    }

    return EXIT_SUCCESS;
}
