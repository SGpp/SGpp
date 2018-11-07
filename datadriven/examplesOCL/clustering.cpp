// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OpFactory.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "sgpp/datadriven/tools/ARFFTools.hpp"

int main() {
  size_t dimension = 10, tiefe = 5, k = 12;  // tiefe 6 for testing
  double lambda = 0.001, treshold = 0.0;
  std::string filename = "../examplesMPI/dataset1_dim8.arff";

  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  start = std::chrono::system_clock::now();

  std::cout << "Loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset data =
    sgpp::datadriven::ARFFTools::readARFFFromFile(filename);
  sgpp::base::DataMatrix& dataset = data.getData();
  dimension = dataset.getNcols();
  std::cout << "Loaded " << dataset.getNcols() << " dimensional dataset with " << dataset.getNrows()
            << " datapoints." << std::endl;

  // Create Grid
  sgpp::base::Grid* grid = sgpp::base::Grid::createLinearGrid(dimension);
  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(tiefe);
  size_t gridsize = grid->getStorage().getSize();
  std::cerr << "Grid created! Number of grid points:     " << gridsize << std::endl;

  sgpp::base::DataVector alpha(gridsize);
  sgpp::base::DataVector result(gridsize);
  alpha.setAll(1.0);

  sgpp::solver::ConjugateGradients* solver = new sgpp::solver::ConjugateGradients(1000, 0.0001);
  sgpp::datadriven::DensityOCLMultiPlatform::OperationDensity* operation_mult =
      sgpp::datadriven::createDensityOCLMultiPlatformConfigured(*grid, dimension, lambda,
                                                                "MyOCLConf.cfg");

  // operation_mult->mult(alpha, result);

  std::cout << "Creating rhs" << std::endl;
  sgpp::base::DataVector b(gridsize);
  operation_mult->generateb(dataset, b);

  // for (size_t i = 0; i < 300; i++)
  //   std::cout << b[i] << " ";
  // std::cout << std::endl;
  // std::ofstream out_rhs("rhs_erg_dim2_depth11.txt");
  // out_rhs.precision(17);
  // for (size_t i = 0; i < gridsize; ++i) {
  //   out_rhs << b[i] << " ";
  // }
  // out_rhs.close();

  std::cout << "Creating alpha" << std::endl;
  solver->solve(*operation_mult, alpha, b, false, true);
  double max = alpha.max();
  double min = alpha.min();
  for (size_t i = 0; i < gridsize; i++) alpha[i] = alpha[i] * 1.0 / (max - min);

  std::ofstream out_alpha("alpha_erg_dim2_depth11.txt");
  out_alpha.precision(17);
  for (size_t i = 0; i < gridsize; ++i) {
    out_alpha << alpha[i] << " ";
  }
  out_alpha.close();
  std::cout << "Starting graph creation..." << std::endl;
  sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL* operation_graph =
      sgpp::datadriven::createNearestNeighborGraphConfigured(dataset, k, dimension,
                                                             "MyOCLConf.cfg");
  std::vector<int> graph(dataset.getNrows() * k);
  operation_graph->create_graph(graph);

  std::ofstream out("graph_erg_dim2_depth11.txt");
  // for (size_t i = 0; i < dataset.getNrows(); ++i) {
  //   for (size_t j = 0; j < k; ++j) {
  //     out << graph[i * k + j] << " ";
  //   }
  //   out << std::endl;
  // }
  // out.close();

  std::cout << "Starting graph pruning" << std::endl;
  sgpp::datadriven::DensityOCLMultiPlatform::OperationPruneGraphOCL* operation_prune =
      sgpp::datadriven::pruneNearestNeighborGraphConfigured(*grid, dimension, alpha, dataset,
                                                            treshold, k, "MyOCLConf.cfg");
  operation_prune->prune_graph(graph);

  // out.open("graph_pruned_erg_dim2_depth11.txt");
  // for (size_t i = 0; i < dataset.getNrows(); ++i) {
  //   for (size_t j = 0; j < k; ++j) {
  //     out << graph[i * k + j] << " ";
  //   }
  //   out << std::endl;
  // }
  // out.close();
  std::vector<size_t> cluster_assignments =
      sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL::find_clusters(graph, k);
  out.open("cluster_erg.txt");
  for (size_t datapoint : cluster_assignments) {
    out << datapoint << " ";
  }
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cout << "finished computation at " << std::ctime(&end_time)
            << "elapsed time: " << elapsed_seconds.count() << "s\n";
  // cleanup
  delete operation_mult;
  delete solver;
  delete operation_graph;
}
