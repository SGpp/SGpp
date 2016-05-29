// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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

#include <iostream>
#include <vector>
#include <string>

#include "sgpp/datadriven/tools/ARFFTools.hpp"

int main(int argc, char *argv[]) {
  // Init MPI enviroment - always has to be done first
  sgpp::datadriven::clusteringmpi::MPIEnviroment::init(argc, argv, true);

  // Loading dataset
  /*std::string filename = "dataset2_dim2.arff";
  std::cout << "Loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset data =
      sgpp::datadriven::ARFFTools::readARFF(filename);
      sgpp::base::DataMatrix& dataset = data.getData();*/

  // Create Grid
  std::unique_ptr<sgpp::base::Grid> grid = sgpp::base::Grid::createLinearGrid(10);
  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(6);
  size_t gridsize = grid->getStorage().getSize();
  std::cerr << "Grid created! Number of grid points:     " << gridsize << std::endl;

  // Create right hand side vector
  /*sgpp::base::DataVector rhs(gridsize);
  sgpp::datadriven::clusteringmpi::OperationRhsMPI rhs_op(*grid, 2, dataset);
  rhs = rhs_op.create_rhs();
  for (auto i = 0; i < 100; ++i) {
    std::cout << rhs[i] << " ";
    }*/

  // Solve for alpha vector via CG solver
  sgpp::base::DataVector alpha(gridsize);
  sgpp::base::DataVector result(gridsize);
  alpha.setAll(1.0);
  sgpp::solver::ConjugateGradients solver(1000, 0.001);
  sgpp::datadriven::clusteringmpi::OperationDensityMPI mult_op(*grid, 0.001);
  //solver.solve(mult_op, alpha, rhs, false, true);
  double max = alpha.max();
  double min = alpha.min();
  /*for (size_t i = 0; i < gridsize; i++)
    alpha[i] = alpha[i]*1.0/(max-min);*/
  mult_op.mult(alpha, result);
  for (auto i = 0; i < 100; ++i) {
    std::cout << result[i] << " ";
  }

  // Create and prune knn graph
  /*sgpp::datadriven::clusteringmpi::OperationCreatePrunedGraph prune_op(*grid, alpha,  dataset, 12);
  std::vector<int> pruned_graph = prune_op.createPrunedGraph(0.7);
  for (size_t i = 0; i < 100; ++i) {
    for (size_t node = 0; node < 12; ++node) {
      std::cout << pruned_graph[i * 12 + node] << " ";
    }
    std::cout << "\n";
    }*/

  // Cleanup MPI enviroment
  sgpp::datadriven::clusteringmpi::MPIEnviroment::release();
  return 0;
}
