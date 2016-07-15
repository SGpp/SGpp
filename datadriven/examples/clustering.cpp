// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OpFactory.hpp>

#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "sgpp/datadriven/tools/ARFFTools.hpp"

int main() {
  size_t dimension = 10, tiefe = 6, k = 12;
  double lambda = 0.000, treshold = 0.7;
  std::string filename = "dataset3_dim10.arff";

  std::cout << "Loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset data =
      sgpp::datadriven::ARFFTools::readARFF(filename);
  sgpp::base::DataMatrix& dataset = data.getData();
  dimension = dataset.getNcols();
  std::cout << "Loaded " << dataset.getNcols() << " dimensional dataset with "
            << dataset.getNrows() << " datapoints." << std::endl;

  /*std::cout << "Size of Grid (3-18): ";
  std::cin >> tiefe;
  std::cout << "Lambda (controlls smoothness of the density function."
            << " 0.01 - 0.0001 recommended.): ";
  std::cin >> lambda;
  std::cout << "k (Number of neighbors for each datapoint. "
            << " 4 - 12 recommended.): ";
  std::cin >> k;
  std::cout << "Treshold (nodes and edges will be removed if their density is below this "
            << "treshold. 0 - 0.2 recommended.): ";
            std::cin >> treshold;*/
  // Create Grid
  std::unique_ptr<sgpp::base::Grid> grid = sgpp::base::Grid::createLinearGrid(dimension);
  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(tiefe);
  size_t gridsize = grid->getStorage().getSize();
  std::cerr << "Grid created! Number of grid points:     " << gridsize << std::endl;

  sgpp::base::DataVector alpha(gridsize);
  sgpp::base::DataVector result(gridsize);
  alpha.setAll(1.0);

  sgpp::solver::ConjugateGradients *solver = new sgpp::solver::ConjugateGradients(1000, 0.001);
  sgpp::datadriven::DensityOCLMultiPlatform::OperationDensityOCL* operation_mult =
      sgpp::datadriven::createDensityOCLMultiPlatformConfigured(*grid, dimension, lambda,
      "MyOCLConf.cfg");

  operation_mult->mult(alpha, result);
  for (auto i = 0; i < 100; ++i) {
    std::cout << result[i] << " ";
  }

  /*std::cout << "Creating rhs" << std::endl;
  sgpp::base::DataVector b(gridsize);
  operation_mult->generateb(dataset, b);
  for (size_t i = 0; i < 300; i++)
    std::cout << b[i] << " ";
  std::cout << std::endl;

  std::cout << "Creating alpha" << std::endl;
  solver->solve(*operation_mult, alpha, b, false, true);
  double max = alpha.max();
  double min = alpha.min();
  for (size_t i = 0; i < gridsize; i++)
  alpha[i] = alpha[i]*1.0/(max-min);

  std::cout << "Starting graph creation..." << std::endl;
  sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL* operation_graph =
      sgpp::datadriven::createNearestNeighborGraphConfigured(dataset, k, dimension,
                                                             "MyOCLConf.cfg");
  std::vector<int> graph(dataset.getNrows()*k);
  operation_graph->create_graph(graph);

  std::ofstream out1("graph_modern1.txt");
  for (size_t i = 0; i < dataset.getNrows(); ++i) {
    for (auto j = 0; j < 12; ++j) {
      out1 << graph[i * 12 + j] << " ";
    }
    out1 << std::endl;
  }
  out1.close();



  std::cout << "Starting graph pruning" << std::endl;
  sgpp::datadriven::DensityOCLMultiPlatform::OperationPruneGraphOCL* operation_prune =
      sgpp::datadriven::pruneNearestNeighborGraphConfigured(*grid, dimension, alpha, dataset,
                                                            treshold, k, "MyOCLConf.cfg");
  operation_prune->prune_graph(graph);

  sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL::find_clusters(graph, k);
  std::ofstream out("graph_modern2.txt");
  for (size_t i = 0; i < dataset.getNrows(); ++i) {
    for (auto j = 0; j < 12; ++j) {
      out << graph[i * 12 + j] << " ";
    }
    out << std::endl;
    }*/
  // cleanup
  delete operation_mult;
  delete solver;
  //delete operation_graph;
}
