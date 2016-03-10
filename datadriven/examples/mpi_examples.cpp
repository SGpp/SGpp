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

#include <iostream>
#include <vector>
#include <string>

#include "sgpp/datadriven/tools/ARFFTools.hpp"

int main(int argc, char *argv[]) {
  // Init MPI enviroment - always has to be done first
  MPIEnviroment::init(argc, argv, true);

  // Loading dataset
  std::string filename = "dataset2_dim2.arff";
  std::cout << "Loading file: " << filename << std::endl;
  SGPP::datadriven::Dataset data =
      SGPP::datadriven::ARFFTools::readARFF(filename);
  SGPP::base::DataMatrix& dataset = data.getData();
  size_t dimensions = dataset.getNcols();

  // Create knn graph
  OperationCreateGraphMPI op(dataset, dimensions, 12);
  std::vector<int> graph = op.create_graph();
  for (size_t i = 0; i < 100; ++i) {
    for (size_t node = 0; node < 12; ++node) {
      std::cout << graph[i * 12 + node] << " ";
    }
    std::cout << "" << "\n";
  }

  // Cleanup MPI enviroment
  MPIEnviroment::release();
  return 0;
}
