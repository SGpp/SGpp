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

int main(int argc, char *argv[]) {

  // Init MPI enviroment - always has to be done first
  sgpp::datadriven::clusteringmpi::MPIEnviroment::init(argc, argv, true);

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
  sgpp::datadriven::clusteringmpi::MPIEnviroment::connect_nodes(testnode);

  // sgpp::datadriven::clusteringmpi::OperationDummy dumdum;
  // dumdum.start_operation();

  // Create Grid
  sgpp::base::Grid *grid = sgpp::base::Grid::createLinearGrid(2);
  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(11);
  size_t gridsize = grid->getStorage().getSize();
  std::cerr << "Grid created! Number of grid points:     " << gridsize << std::endl;

  sgpp::datadriven::clusteringmpi::OperationDensityMultMPI mult_op(*grid, 0.001);

  sgpp::base::DataVector alpha(gridsize);
  sgpp::base::DataVector result(gridsize);
  alpha.setAll(1.0);
  mult_op.mult(alpha, result);
  std::cout << std::endl << std::endl;
  for (size_t i = 0; i < 100; ++i) {
    std::cout << result[i] << " ";
  }
  std::cin.get();
  std::cout << std::endl << std::endl;


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

  // Create and prune knn graph
  std::cin.get();
  sgpp::datadriven::clusteringmpi::OperationPrunedGraphCreationMPI graph_op(*grid, alpha,
                                                                            dataset, 12, 0.7);
  std::cin.get();
  std::vector<int> knn_graph;
  graph_op.create_graph(knn_graph);
  for (size_t i = 0; i < 100; ++i) {
    for (size_t node = 0; node < 12; ++node) {
      std::cout << knn_graph[i * 12 + node] << " ";
    }
    std::cout << "\n";
    }

  // Cleanup MPI enviroment
  sgpp::datadriven::clusteringmpi::MPIEnviroment::release();
  return 0;
}
