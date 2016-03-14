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

#include <iostream>
#include <vector>
#include <string>

#include "sgpp/datadriven/tools/ARFFTools.hpp"

int main(int argc, char *argv[]) {
  // Init MPI enviroment - always has to be done first
  sgpp::datadriven::clusteringmpi::MPIEnviroment::init(argc, argv, true);

  // Loading dataset
  std::string filename = "dataset2_dim2.arff";
  std::cout << "Loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset data =
      sgpp::datadriven::ARFFTools::readARFF(filename);
  sgpp::base::DataMatrix& dataset = data.getData();
  size_t dimensions = dataset.getNcols();

  // Create Grid
  std::unique_ptr<sgpp::base::Grid> grid = sgpp::base::Grid::createLinearGrid(2);
  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(10);
  size_t gridsize = grid->getStorage().getSize();
  std::cerr << "Grid created! Number of grid points:     " << gridsize << std::endl;

  sgpp::base::DataVector result(gridsize);
  sgpp::datadriven::clusteringmpi::OperationRhsMPI rhs_op(*grid, 2, dataset);
  result = rhs_op.create_rhs();
  for (auto i = 0; i < 100; ++i) {
    std::cout << result[i] << " ";
  }

  // Create knn graph
  sgpp::datadriven::clusteringmpi::OperationCreateGraphMPI op(dataset, dimensions, 12);
  std::vector<int> graph = op.create_graph();
  for (size_t i = 0; i < 100; ++i) {
    for (size_t node = 0; node < 12; ++node) {
      std::cout << graph[i * 12 + node] << " ";
    }
    std::cout << "" << "\n";
  }

  // Cleanup MPI enviroment
  sgpp::datadriven::clusteringmpi::MPIEnviroment::release();
  return 0;
}
