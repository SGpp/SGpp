// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#include <unistd.h>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/MPIEnviroment.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationCreateGraphMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationRhsMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationDensityMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationCreatePrunedGraphMPI.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <sgpp/datadriven/operation/hash/OperationMPI/OperationGridBaseMPI.hpp>

#include <iostream>
#include <vector>
#include <string>

#include "sgpp/datadriven/tools/ARFFTools.hpp"

// arg 1: grid dimensions
// arg 2: grid level
// arg 3: packagesize - master
// arg 4: packagesize - leutant
// arg 5: output file
int main(int argc, char *argv[]) {
  // Init MPI enviroment - always has to be done first
  sgpp::datadriven::clusteringmpi::MPIEnviroment::init(argc, argv, true);

  if (argc != 6) {
    std::cout << "Wrong arguments - consult source file for more informations!" << std::endl;
    sgpp::datadriven::clusteringmpi::MPIEnviroment::release();
    return 0;
  }

  int dim = std::stoi(argv[1]);
  int level = std::stoi(argv[2]);
  int packagesize_master = std::stoi(argv[3]);
  int packagesize_leut = std::stoi(argv[4]);

  // Create Sample config
  sgpp::base::OperationConfiguration conf = sgpp::datadriven::clusteringmpi::
                                            MPIEnviroment::createMPIConfiguration(2, 4);
  conf.serialize("testconf.cfg");
  sgpp::base::OCLOperationConfiguration ocl_conf("MyOCLConf.cfg");
  sgpp::base::OperationConfiguration conf_ocl = sgpp::datadriven::clusteringmpi::
                                                MPIEnviroment::createMPIConfiguration(2,
                                                                                      ocl_conf);
  conf_ocl.serialize("testconf2.cfg");

  // MPI_Init(&argc, &argv);
  sgpp::base::OperationConfiguration testnode("MPIConf2.cfg");
  testnode["PREFERED_PACKAGESIZE"].setInt(packagesize_master);
  sgpp::datadriven::clusteringmpi::MPIEnviroment::connect_nodes(testnode);

  // sgpp::datadriven::clusteringmpi::OperationDummy dumdum;
  // dumdum.start_operation();

  // Create Grid
  sgpp::base::Grid *grid = sgpp::base::Grid::createLinearGrid(dim);
  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(level);
  size_t gridsize = grid->getStorage().getSize();
  std::cerr << "Grid created! Number of grid points:     " << gridsize << std::endl;
  std::cin.get();

  sgpp::datadriven::clusteringmpi::OperationDensityMultMPI mult_op(*grid, 0.001);

  sgpp::base::DataVector alpha(gridsize);
  sgpp::base::DataVector result(gridsize);
  alpha.setAll(1.0);
  mult_op.mult(alpha, result);
  std::cout << std::endl << std::endl;
  for (size_t i = 0; i < 100; ++i) {
    std::cout << result[i] << " ";
  }
  std::cout << std::endl << std::endl;

  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  start = std::chrono::system_clock::now();
  //sgpp::datadriven::clusteringmpi::OperationGridMethod test(testnode, *grid, "grid_dummy");
  // Loading dataset
  std::string filename = "dataset2_dim2.arff";
  std::cout << "Loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset data =
      sgpp::datadriven::ARFFTools::readARFF(filename);
  sgpp::base::DataMatrix& dataset = data.getData();

  // Create right hand side vector
  sgpp::base::DataVector rhs(gridsize);
  sgpp::datadriven::clusteringmpi::OperationDensityRhsMPI rhs_op(*grid, dataset);
  rhs_op.generate_b(rhs);
  for (auto i = 0; i < 100; ++i) {
    std::cout << rhs[i] << " ";
  }

  // Solve for alpha vector via CG solver
  alpha.setAll(1.0);
  sgpp::solver::ConjugateGradients solver(1000, 0.001);
  //sgpp::datadriven::clusteringmpi::OperationDensityMultMPI mult_op(*grid, 0.001);
  solver.solve(mult_op, alpha, rhs, false, true);
  double max = alpha.max();
  double min = alpha.min();
  for (size_t i = 0; i < gridsize; i++)
    alpha[i] = alpha[i]*1.0/(max-min);
  end = std::chrono::system_clock::now();

  // Calc time
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cout << "finished computation at " << std::ctime(&end_time)
            << "elapsed time: " << elapsed_seconds.count() << "s\n";
  // Write result into file
  std::ofstream ofs;
  ofs.open (argv[5], std::ofstream::out | std::ofstream::app);
  ofs << dim << ";" << level << ";" << packagesize_master << ";" << packagesize_leut
      << ";" << elapsed_seconds.count() << ";"
      << gridsize % ((sgpp::datadriven::clusteringmpi::MPIEnviroment::get_node_count() -1) *
                     packagesize_master)<< std::endl;
  ofs.close();

  // Create and prune knn graph
  /*sgpp::datadriven::clusteringmpi::OperationPrunedGraphCreationMPI graph_op(*grid, alpha,
    dataset, 12, 0.7);
    std::vector<int> knn_graph;
    graph_op.create_graph(knn_graph);
    std::cout << "knn graph size: " << knn_graph.size() / 12 << std::endl;
    for (size_t i = 0; i < 100; ++i) {
    for (size_t node = 0; node < 12; ++node) {
    std::cout << knn_graph[i * 12 + node] << " ";
    }
    std::cout << "\n";
    }
    sgpp::datadriven::DensityOCLMultiPlatform::
    OperationCreateGraphOCL::find_clusters(knn_graph, 12);
  */

  // Cleanup MPI enviroment
  sgpp::datadriven::clusteringmpi::MPIEnviroment::release();
  return 0;
}
