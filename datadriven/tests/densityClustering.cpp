// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#if USE_OCL == 1

#define BOOST_TEST_DYN_LINK
#include <zlib.h>
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <vector>
#include "sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OperationCreateGraphOCLSingleDevice.hpp"
#include "sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp"
#include "sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OpFactory.hpp"
#include "sgpp/solver/sle/ConjugateGradients.hpp"

#include "sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/base/tools/ConfigurationParameters.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/globaldef.hpp"
#include "test_datadrivenCommon.hpp"
void multiply_and_test(sgpp::base::OCLOperationConfiguration *parameters,
                       std::vector<double> &mult_optimal_result,
                       std::shared_ptr<sgpp::base::OCLManagerMultiPlatform> manager,
                       sgpp::base::Grid &grid) {
  size_t gridsize = grid.getStorage().getSize();
  // Create vectors for multiplication
  sgpp::base::DataVector alpha(gridsize);
  sgpp::base::DataVector result(gridsize);
  alpha.setAll(1.0);
  // Create operation
  sgpp::datadriven::DensityOCLMultiPlatform::OperationDensity *mult_operation =
      new sgpp::datadriven::DensityOCLMultiPlatform::OperationDensityOCLMultiPlatform<double>(
          grid, 2, manager, parameters, 0.001, 0, 0);
  // Execute multiplication
  mult_operation->mult(alpha, result);
  // Compare results with optimal results
  for (size_t i = 0; i < gridsize; ++i) {
    BOOST_CHECK_CLOSE(mult_optimal_result[i], result[i], 0.001);
  }
  delete mult_operation;
}

BOOST_AUTO_TEST_SUITE(TestClusteringOpenCL)

BOOST_AUTO_TEST_CASE(DensityMultiplicationOpenCL) {
  // Load correct results for comparison
  std::vector<double> mult_optimal_result;
  std::ifstream mult_in("datadriven/datasets/clustering_test_data/mult_erg_dim2_depth11.txt");
  if (mult_in) {
    double value;
    while (mult_in >> value) mult_optimal_result.push_back(value);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("Density multiplication result file is missing!"));
  }
  mult_in.close();

  // Create OCL configuration
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();
  sgpp::datadriven::DensityOCLMultiPlatform::OperationDensity::load_default_parameters(
      parameters.get());

  // Create OpenCL Manager
  std::shared_ptr<sgpp::base::OCLManagerMultiPlatform> manager;
  manager = std::make_shared<sgpp::base::OCLManagerMultiPlatform>(false);

  // Create grid for test scenario
  sgpp::base::Grid *grid = sgpp::base::Grid::createLinearGrid(2);
  sgpp::base::GridGenerator &gridGen = grid->getGenerator();
  gridGen.regular(11);

  std::cout << "Testing default kernel configuration..." << std::endl;
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Testing with preprocessed positions..." << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "multdensity";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("PREPROCESS_POSITIONS", true);
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Testing with preprocessed positions and ignored flags..." << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "multdensity";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("PREPROCESS_POSITIONS", true);
      kernelNode.replaceIDAttr("USE_FABS", true);
      kernelNode.replaceIDAttr("USE_IMPLICIT", true);
      kernelNode.replaceIDAttr("USE_LEVEL_CACHE", true);
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "multdensity";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("PREPROCESS_POSITIONS", true);
      kernelNode.replaceIDAttr("USE_FABS", true);
      kernelNode.replaceIDAttr("USE_IMPLICIT", true);
      kernelNode.replaceIDAttr("USE_LESS_OPERATIONS", false);
      kernelNode.replaceIDAttr("USE_LEVEL_CACHE", false);
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Testing with branchless mutliplication kernel..." << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "multdensity";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("PREPROCESS_POSITIONS", false);
      kernelNode.replaceIDAttr("USE_FABS", false);
      kernelNode.replaceIDAttr("USE_IMPLICIT", false);
      kernelNode.replaceIDAttr("USE_LESS_OPERATIONS", false);
      kernelNode.replaceIDAttr("USE_LEVEL_CACHE", false);
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Testing with branchless mutliplication kernel with fabs modifications..."
            << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "multdensity";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("PREPROCESS_POSITIONS", false);
      kernelNode.replaceIDAttr("USE_FABS", true);
      kernelNode.replaceIDAttr("USE_IMPLICIT", false);
      kernelNode.replaceIDAttr("USE_LESS_OPERATIONS", false);
      kernelNode.replaceIDAttr("USE_LEVEL_CACHE", false);
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);
  std::cout << "Testing with branchless mutliplication kernel with level cache..." << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "multdensity";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("PREPROCESS_POSITIONS", false);
      kernelNode.replaceIDAttr("USE_FABS", true);
      kernelNode.replaceIDAttr("USE_IMPLICIT", false);
      kernelNode.replaceIDAttr("USE_LESS_OPERATIONS", false);
      kernelNode.replaceIDAttr("USE_LEVEL_CACHE", false);
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Testing with branchless mutliplication kernel with implicit modifications..."
            << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "multdensity";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("PREPROCESS_POSITIONS", false);
      kernelNode.replaceIDAttr("USE_FABS", false);
      kernelNode.replaceIDAttr("USE_IMPLICIT", true);
      kernelNode.replaceIDAttr("USE_LESS_OPERATIONS", false);
      kernelNode.replaceIDAttr("USE_LEVEL_CACHE", false);
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Testing with branchless mutliplication kernel with "
            << "all modifications..." << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "multdensity";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("PREPROCESS_POSITIONS", false);
      kernelNode.replaceIDAttr("USE_FABS", true);
      kernelNode.replaceIDAttr("USE_IMPLICIT", true);
      kernelNode.replaceIDAttr("USE_LESS_OPERATIONS", false);
      kernelNode.replaceIDAttr("USE_LEVEL_CACHE", true);
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Testing default multiplcation kernel with fabs modifications..." << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "multdensity";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("PREPROCESS_POSITIONS", false);
      kernelNode.replaceIDAttr("USE_FABS", true);
      kernelNode.replaceIDAttr("USE_IMPLICIT", false);
      kernelNode.replaceIDAttr("USE_LESS_OPERATIONS", true);
      kernelNode.replaceIDAttr("USE_LEVEL_CACHE", false);
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Testing default multiplcation kernel with implicit modifications..." << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "multdensity";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("PREPROCESS_POSITIONS", false);
      kernelNode.replaceIDAttr("USE_FABS", false);
      kernelNode.replaceIDAttr("USE_IMPLICIT", true);
      kernelNode.replaceIDAttr("USE_LESS_OPERATIONS", true);
      kernelNode.replaceIDAttr("USE_LEVEL_CACHE", false);
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Testing default multiplcation kernel with level cache..." << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "multdensity";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("PREPROCESS_POSITIONS", false);
      kernelNode.replaceIDAttr("USE_FABS", false);
      kernelNode.replaceIDAttr("USE_IMPLICIT", false);
      kernelNode.replaceIDAttr("USE_LESS_OPERATIONS", true);
      kernelNode.replaceIDAttr("USE_LEVEL_CACHE", true);
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Testing default multiplcation kernel with all modifications..." << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "multdensity";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("PREPROCESS_POSITIONS", false);
      kernelNode.replaceIDAttr("USE_FABS", true);
      kernelNode.replaceIDAttr("USE_IMPLICIT", true);
      kernelNode.replaceIDAttr("USE_LESS_OPERATIONS", true);
      kernelNode.replaceIDAttr("USE_LEVEL_CACHE", true);
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);
}

BOOST_AUTO_TEST_CASE(DensityAlphaSolver) {
  // Load correct results for comparison
  std::vector<double> alpha_optimal_result;
  std::ifstream alpha_in("datadriven/datasets/clustering_test_data/alpha_erg_dim2_depth11.txt");
  if (alpha_in) {
    double value;
    while (alpha_in >> value) alpha_optimal_result.push_back(value);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("Alpha result file is missing!"));
  }
  alpha_in.close();

  // Create grid for test scenario
  sgpp::base::Grid *grid = sgpp::base::Grid::createLinearGrid(2);
  sgpp::base::GridGenerator &gridGen = grid->getGenerator();
  gridGen.regular(11);
  size_t gridsize = grid->getStorage().getSize();

  // Load rhs vector
  sgpp::base::DataVector b(gridsize);
  std::ifstream rhs_in("datadriven/datasets/clustering_test_data/rhs_erg_dim2_depth11.txt");
  if (rhs_in) {
    double value;
    int counter = 0;
    while (rhs_in >> value) {
      b[counter] = value;
      counter++;
    }
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("Density rhs result file is missing!"));
  }
  rhs_in.close();

  // Create OCL configuration
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();
  sgpp::datadriven::DensityOCLMultiPlatform::OperationDensity::load_default_parameters(
      parameters.get());

  // Create OpenCL Manager
  std::shared_ptr<sgpp::base::OCLManagerMultiPlatform> manager;
  manager = std::make_shared<sgpp::base::OCLManagerMultiPlatform>(false);

  // Create operation
  sgpp::datadriven::DensityOCLMultiPlatform::OperationDensity *mult_operation =
      new sgpp::datadriven::DensityOCLMultiPlatform::OperationDensityOCLMultiPlatform<double>(
          *grid, 2, manager, parameters.get(), 0.001, 0, 0);

  // Create solver
  sgpp::solver::ConjugateGradients *solver = new sgpp::solver::ConjugateGradients(100, 0.001);
  sgpp::base::DataVector alpha(gridsize);
  alpha.setAll(1.0);

  // Solve
  solver->solve(*mult_operation, alpha, b, false, false);
  // Scaling
  double max = alpha.max();
  double min = alpha.min();
  for (size_t i = 0; i < gridsize; i++) alpha[i] = alpha[i] * 1.0 / (max - min);

  // Compare results with correct results
  for (size_t i = 0; i < gridsize; ++i) {
    BOOST_CHECK_CLOSE(alpha_optimal_result[i], alpha[i], 1.0);
  }

  // Cleanup
  delete mult_operation;
  delete solver;
}

BOOST_AUTO_TEST_CASE(DensityRHSOpenCL) {
  // Load correct results for comparison
  std::vector<double> rhs_optimal_result;
  std::ifstream rhs_in("datadriven/datasets/clustering_test_data/rhs_erg_dim2_depth11.txt");
  if (rhs_in) {
    double value;
    while (rhs_in >> value) rhs_optimal_result.push_back(value);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("Density rhs result file is missing!"));
  }
  rhs_in.close();

  // Create OCL configuration
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();
  sgpp::datadriven::DensityOCLMultiPlatform::OperationDensity::load_default_parameters(
      parameters.get());

  // Create OpenCL Manager
  std::shared_ptr<sgpp::base::OCLManagerMultiPlatform> manager;
  manager = std::make_shared<sgpp::base::OCLManagerMultiPlatform>(false);

  // Create grid for test scenario
  sgpp::base::Grid *grid = sgpp::base::Grid::createLinearGrid(2);
  sgpp::base::GridGenerator &gridGen = grid->getGenerator();
  gridGen.regular(11);
  size_t gridsize = grid->getStorage().getSize();

  // Load dataset for test scenario
  sgpp::datadriven::Dataset data =
    sgpp::datadriven::ARFFTools::readARFFFromFile(
      "datadriven/datasets/clustering_test_data/clustering_testdataset_dim2.arff");
  sgpp::base::DataMatrix &dataset = data.getData();

  // Create operation
  sgpp::datadriven::DensityOCLMultiPlatform::OperationDensity *operation_rhs =
      new sgpp::datadriven::DensityOCLMultiPlatform::OperationDensityOCLMultiPlatform<double>(
          *grid, 2, manager, parameters.get(), 0.001, 0, 0);

  std::cout << "Testing rhs kernel ..." << std::endl;
  sgpp::base::DataVector b(gridsize);
  operation_rhs->generateb(dataset, b);
  for (size_t i = 0; i < gridsize; ++i) {
    BOOST_CHECK_CLOSE(rhs_optimal_result[i], b[i], 0.001);
  }
  delete operation_rhs;
}

BOOST_AUTO_TEST_CASE(KNNGraphOpenCL) {
  // Load correct results for comparison
  std::vector<int> graph_optimal_result;
  std::ifstream graph_in("datadriven/datasets/clustering_test_data/graph_erg_dim2_depth11.txt");
  if (graph_in) {
    int value;
    while (graph_in >> value) graph_optimal_result.push_back(value);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("knn graph result file is missing!"));
  }
  graph_in.close();

  // Create OCL configuration
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();
  sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL::load_default_parameters(
      parameters.get());

  // Create OpenCL Manager
  std::shared_ptr<sgpp::base::OCLManagerMultiPlatform> manager;
  manager = std::make_shared<sgpp::base::OCLManagerMultiPlatform>(false);

  // Load dataset for test scenario
  sgpp::datadriven::Dataset data =
    sgpp::datadriven::ARFFTools::readARFFFromFile(
      "datadriven/datasets/clustering_test_data/clustering_testdataset_dim2.arff");
  sgpp::base::DataMatrix &dataset = data.getData();

  // Create operation
  size_t k = 8;
  sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL *operation_graph =
      new sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCLSingleDevice<double>(
          dataset, 2, manager, parameters.get(), k, 0, 0);
  // Test graph kernel
  std::vector<int> graph(dataset.getNrows() * k);
  operation_graph->create_graph(graph);
  for (size_t i = 0; i < dataset.getNrows() * k; ++i) {
    BOOST_CHECK(graph_optimal_result[i] == graph[i]);
  }
  delete operation_graph;
  operation_graph = NULL;

  std::cout << "Testing default knn graph kernel with select statements..." << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "connectNeighbors";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("USE_SELECT", true);
    }
  }
  operation_graph =
      new sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCLSingleDevice<double>(
          dataset, 2, manager, parameters.get(), k, 0, 0);
  // Test graph kernel
  operation_graph->create_graph(graph);
  for (size_t i = 0; i < dataset.getNrows() * k; ++i) {
    BOOST_CHECK(graph_optimal_result[i] == graph[i]);
  }

  delete operation_graph;
  operation_graph = NULL;

  std::cout << "Testing default knn graph kernel with local memory..." << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "connectNeighbors";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("USE_SELECT", false);
      kernelNode.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
    }
  }
  operation_graph =
      new sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCLSingleDevice<double>(
          dataset, 2, manager, parameters.get(), k, 0, 0);
  // Test graph kernel
  operation_graph->create_graph(graph);
  for (size_t i = 0; i < dataset.getNrows() * k; ++i) {
    BOOST_CHECK(graph_optimal_result[i] == graph[i]);
  }
  delete operation_graph;
  operation_graph = NULL;

  std::cout << "Testing default knn graph kernel with local memory and select statements..."
            << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "connectNeighbors";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("USE_SELECT", false);
      kernelNode.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
    }
  }
  operation_graph =
      new sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCLSingleDevice<double>(
          dataset, 2, manager, parameters.get(), k, 0, 0);
  // Test graph kernel
  operation_graph->create_graph(graph);
  for (size_t i = 0; i < dataset.getNrows() * k; ++i) {
    BOOST_CHECK(graph_optimal_result[i] == graph[i]);
  }
  delete operation_graph;
}

BOOST_AUTO_TEST_CASE(KNNPruneGraphOpenCL) {
  // Load input
  std::vector<int> graph;
  std::ifstream graph_in("datadriven/datasets/clustering_test_data/graph_erg_dim2_depth11.txt");
  if (graph_in) {
    int value;
    while (graph_in >> value) graph.push_back(value);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("knn graph result file is missing!"));
  }
  graph_in.close();

  // Load optimal results for comparison
  std::vector<int> graph_optimal_result;
  std::ifstream graph_result(
      "datadriven/datasets/clustering_test_data/graph_pruned_erg_dim2_depth11.txt");
  if (graph_result) {
    int value;
    while (graph_result >> value) graph_optimal_result.push_back(value);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("knn pruned graph result file is missing!"));
  }
  graph_result.close();

  // Create grid for test scenario
  sgpp::base::Grid *grid = sgpp::base::Grid::createLinearGrid(2);
  sgpp::base::GridGenerator &gridGen = grid->getGenerator();
  gridGen.regular(11);
  size_t gridsize = grid->getStorage().getSize();

  // Load optimal results for comparison
  sgpp::base::DataVector alpha(gridsize);
  std::ifstream alpha_result(
      "datadriven/datasets/clustering_test_data/alpha_erg_dim2_depth11.txt");
  if (alpha_result) {
    size_t counter = 0;
    double value;
    while (alpha_result >> value) {
      alpha[counter] = value;
      counter++;
    }
    BOOST_CHECK(counter == gridsize);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("knn pruned graph result file is missing!"));
  }
  alpha_result.close();

  // Create OCL configuration
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();
  sgpp::datadriven::DensityOCLMultiPlatform::OperationPruneGraphOCL::load_default_parameters(
      parameters.get());

  // Create OpenCL Manager
  std::shared_ptr<sgpp::base::OCLManagerMultiPlatform> manager;
  manager = std::make_shared<sgpp::base::OCLManagerMultiPlatform>(false);

  // Load dataset for test scenario
  sgpp::datadriven::Dataset data =
    sgpp::datadriven::ARFFTools::readARFFFromFile(
      "datadriven/datasets/clustering_test_data/clustering_testdataset_dim2.arff");
  sgpp::base::DataMatrix &dataset = data.getData();

  // Create operation
  sgpp::datadriven::DensityOCLMultiPlatform::OperationPruneGraphOCL *operation_prune =
      new sgpp::datadriven::DensityOCLMultiPlatform::OperationPruneGraphOCLMultiPlatform<double>(
          *grid, alpha, dataset, 2, manager, parameters.get(), 0.2, 8, 0, 0);

  std::cout << "Testing knn prune kernel ..." << std::endl;
  operation_prune->prune_graph(graph);
  for (size_t i = 0; i < gridsize; ++i) {
    BOOST_CHECK(graph[i] == graph_optimal_result[i]);
  }
  delete operation_prune;
}

BOOST_AUTO_TEST_CASE(KNNClusterSearch) {
  // Load input
  std::vector<int> graph;
  std::ifstream graph_in(
      "datadriven/datasets/clustering_test_data/graph_pruned_erg_dim2_depth11.txt");
  if (graph_in) {
    int value;
    while (graph_in >> value) graph.push_back(value);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("Pruned knn graph result file is missing!"));
  }
  graph_in.close();

  std::vector<size_t> optimal_cluster_assignement;
  std::ifstream assignement_in("datadriven/datasets/clustering_test_data/cluster_erg.txt");
  if (assignement_in) {
    size_t value;
    while (assignement_in >> value) optimal_cluster_assignement.push_back(value);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("Pruned knn graph result file is missing!"));
  }
  assignement_in.close();

  std::vector<size_t> cluster_assignement =
      sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL::find_clusters(graph, 8);
  BOOST_CHECK(optimal_cluster_assignement.size() == cluster_assignement.size());
  if (optimal_cluster_assignement.size() == cluster_assignement.size()) {
    for (size_t i = 0; i < cluster_assignement.size(); ++i) {
      BOOST_CHECK(optimal_cluster_assignement[i] == cluster_assignement[i]);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
#endif
