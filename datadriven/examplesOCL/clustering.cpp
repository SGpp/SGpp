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
  size_t dimension = 10, tiefe = 6, k = 8;
  double lambda = 0.0001, treshold = 0.2;
  std::string filename = "clustering_testdataset_dim2.arff";

  std::cout << "Loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset data =
      sgpp::datadriven::ARFFTools::readARFF(filename);
  sgpp::base::DataMatrix& dataset = data.getData();
  dimension = dataset.getNcols();
  std::cout << "Loaded " << dataset.getNcols() << " dimensional dataset with "
            << dataset.getNrows() << " datapoints." << std::endl;

  // Create Grid
  sgpp::base::Grid *grid = sgpp::base::Grid::createLinearGrid(2);
  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(11);
  size_t gridsize = grid->getStorage().getSize();
  std::cerr << "Grid created! Number of grid points:     " << gridsize << std::endl;

  sgpp::base::DataVector alpha(gridsize);
  sgpp::base::DataVector result(gridsize);
  alpha.setAll(1.0);

  sgpp::solver::ConjugateGradients *solver = new sgpp::solver::ConjugateGradients(1000, 0.001);
  sgpp::datadriven::DensityOCLMultiPlatform::OperationDensityOCL* operation_mult =
      sgpp::datadriven::createDensityOCLMultiPlatformConfigured(*grid, dimension, 0.001,
                                                                "MyOCLConf.cfg");

  operation_mult->mult(alpha, result);

  std::ofstream out_mult("mult_erg_dim2_depth11.txt");
  out_mult.precision(17);
  for (size_t i = 0; i < gridsize; ++i) {
    out_mult << result[i] << " ";
  }
  out_mult.close();

  std::cout << "Creating rhs" << std::endl;
  sgpp::base::DataVector b(gridsize);
  operation_mult->generateb(dataset, b);
  for (size_t i = 0; i < 300; i++)
    std::cout << b[i] << " ";
  std::cout << std::endl;
  std::ofstream out_rhs("rhs_erg_dim2_depth11.txt");
  out_rhs.precision(17);
  for (size_t i = 0; i < gridsize; ++i) {
    out_rhs << b[i] << " ";
  }
  out_rhs.close();

  std::cout << "Creating alpha" << std::endl;
  solver->solve(*operation_mult, alpha, b, false, true);
  double max = alpha.max();
  double min = alpha.min();
  for (size_t i = 0; i < gridsize; i++)
    alpha[i] = alpha[i]*1.0/(max-min);

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
  std::vector<int> graph(dataset.getNrows()*k);
  operation_graph->create_graph(graph);

  std::ofstream out("graph_erg_dim2_depth11.txt");
  for (size_t i = 0; i < dataset.getNrows(); ++i) {
    for (size_t j = 0; j < k; ++j) {
      out << graph[i * k + j] << " ";
    }
    out << std::endl;
  }
  out.close();

  std::cout << "Starting graph pruning" << std::endl;
  sgpp::datadriven::DensityOCLMultiPlatform::OperationPruneGraphOCL* operation_prune =
      sgpp::datadriven::pruneNearestNeighborGraphConfigured(*grid, dimension, alpha, dataset,
                                                            treshold, k, "MyOCLConf.cfg");
  operation_prune->prune_graph(graph);

  out.open("graph_pruned_erg_dim2_depth11.txt");
  for (size_t i = 0; i < dataset.getNrows(); ++i) {
    for (size_t j = 0; j < k; ++j) {
      out << graph[i * k + j] << " ";
    }
    out << std::endl;
  }
  out.close();
  std::vector<size_t> cluster_assignments =
      sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL::find_clusters(graph, k);
  out.open("cluster_erg.txt");
  for (size_t datapoint : cluster_assignments) {
    out << datapoint << " ";
  }
  // cleanup
  delete operation_mult;
  delete solver;
  delete operation_graph;
}
