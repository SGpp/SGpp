// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#if USE_OCL == 1

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <zlib.h>

#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OperationCreateGraphOCLSingleDevice.hpp>
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OpFactory.hpp>
#include <random>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "test_datadrivenCommon.hpp"
#include "sgpp/globaldef.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "sgpp/base/tools/ConfigurationParameters.hpp"
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
  sgpp::datadriven::DensityOCLMultiPlatform::OperationDensityOCL* mult_operation =
      new sgpp::datadriven::DensityOCLMultiPlatform::
      OperationDensityOCLMultiPlatform<double>(grid, 2, manager, parameters, 0.001,
                                               0, 0);
  // Execute multiplication
  mult_operation->mult(alpha, result);
  // Compare results with optimal results
  for (size_t i = 0; i < gridsize; ++i) {
    BOOST_CHECK_CLOSE(mult_optimal_result[i], result[i], 0.001);
  }
  delete mult_operation;
}

BOOST_AUTO_TEST_SUITE(TestClusteringOpenCL)

BOOST_AUTO_TEST_CASE(DensityMultiplicationOpenCL)  {
  // Load correct results for comparison
  std::vector<double> mult_optimal_result;
  std::ifstream mult_in( "/media/DatenPartition/Programmieren/ipvs/SGpp/datadriven/tests/datasets/mult_erg_dim2_depth11.txt");
  if (mult_in) {
    double value;
    while (mult_in >> value)
      mult_optimal_result.push_back(value);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("Density multiplication result file is missing!"));
  }

  // Create OCL configuration
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();
  sgpp::datadriven::DensityOCLMultiPlatform::
      OperationDensityOCL::load_default_parameters(parameters.get());

  // Create OpenCL Manager
  std::shared_ptr<sgpp::base::OCLManagerMultiPlatform> manager;
  manager = std::make_shared<sgpp::base::OCLManagerMultiPlatform>(true);

  // Create grid for test scenario
  sgpp::base::Grid *grid = sgpp::base::Grid::createLinearGrid(2);
  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(11);

  std::cout << "Test default kernel configuration..." << std::endl;
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Test with preprocessed positions:" << std::endl;
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

  std::cout << "Test with preprocessed positions and ignored flags..." << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "multdensity";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("PREPROCESS_POSITIONS", true);
      kernelNode.replaceIDAttr("USE_FABS", true);
      kernelNode.replaceIDAttr("USE_IMPLICIT", true);
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Test with branchless mutliplication kernel..." << std::endl;
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
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Test default multiplcation kernel with fabs modifications..." << std::endl;
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
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Test default multiplcation kernel with implicit modifications..." << std::endl;
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
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Test default multiplcation kernel with implicit and fabs modifications..."
            << std::endl;
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
    }
  }
  multiply_and_test(parameters.get(), mult_optimal_result, manager, *grid);

  std::cout << "Test default multiplcation kernel with level cache..."
            << std::endl;
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
}

BOOST_AUTO_TEST_CASE(DensityRHSOpenCL)  {
  // Load correct results for comparison
  std::vector<double> rhs_optimal_result;
  std::ifstream rhs_in( "/media/DatenPartition/Programmieren/ipvs/SGpp/datadriven/tests/datasets/rhs_erg_dim2_depth11.txt");
  if (rhs_in) {
    double value;
    while (rhs_in >> value)
      rhs_optimal_result.push_back(value);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("Density rhs result file is missing!"));
  }

  // Create OCL configuration
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();
  sgpp::datadriven::DensityOCLMultiPlatform::
      OperationDensityOCL::load_default_parameters(parameters.get());

  // Create OpenCL Manager
  std::shared_ptr<sgpp::base::OCLManagerMultiPlatform> manager;
  manager = std::make_shared<sgpp::base::OCLManagerMultiPlatform>(true);

  // Create grid for test scenario
  sgpp::base::Grid *grid = sgpp::base::Grid::createLinearGrid(2);
  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(11);
  size_t gridsize = grid->getStorage().getSize();

  // Load dataset for test scenario
  sgpp::datadriven::Dataset data =
      sgpp::datadriven::ARFFTools::
      readARFF("/media/DatenPartition/Programmieren/ipvs/SGpp/datadriven/tests/datasets/clustering_testdataset_dim2.arff");
  sgpp::base::DataMatrix& dataset = data.getData();

  // Create operation
  sgpp::datadriven::DensityOCLMultiPlatform::OperationDensityOCL* operation_rhs =
      new sgpp::datadriven::DensityOCLMultiPlatform::
      OperationDensityOCLMultiPlatform<double>(*grid, 2, manager, parameters.get(), 0.001,
                                               0, 0);

  std::cout << "Test rhs kernel ..." << std::endl;
  sgpp::base::DataVector b(gridsize);
  operation_rhs->generateb(dataset, b);
  for (size_t i = 0; i < gridsize; ++i) {
    BOOST_CHECK_CLOSE(rhs_optimal_result[i], b[i], 0.001);
  }
  delete operation_rhs;
}

BOOST_AUTO_TEST_CASE(KNNGraphOpenCL)  {
  // Load correct results for comparison
  std::vector<int> graph_optimal_result;
  std::ifstream graph_in( "/media/DatenPartition/Programmieren/ipvs/SGpp/datadriven/tests/datasets/graph_erg_dim2_depth11.txt");
  if (graph_in) {
    double value;
    while (graph_in >> value)
      graph_optimal_result.push_back(value);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("knn graph result file is missing!"));
  }

  // Create OCL configuration
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();
  sgpp::datadriven::DensityOCLMultiPlatform::
      OperationCreateGraphOCL::load_default_parameters(parameters.get());

  // Create OpenCL Manager
  std::shared_ptr<sgpp::base::OCLManagerMultiPlatform> manager;
  manager = std::make_shared<sgpp::base::OCLManagerMultiPlatform>(true);

  // Load dataset for test scenario
  sgpp::datadriven::Dataset data =
      sgpp::datadriven::ARFFTools::
      readARFF("/media/DatenPartition/Programmieren/ipvs/SGpp/datadriven/tests/datasets/clustering_testdataset_dim2.arff");
  sgpp::base::DataMatrix& dataset = data.getData();

  // Create operation
  size_t k = 8;
  sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL* operation_graph =
      new sgpp::datadriven::DensityOCLMultiPlatform::
        OperationCreateGraphOCLSingleDevice<double>(dataset, 2, manager,
                                                    parameters.get(), k, 0, 0);
  // Test graph kernel
  std::vector<int> graph(dataset.getNrows() * k);
  operation_graph->create_graph(graph);
  for (size_t i = 0; i < dataset.getNrows() * k; ++i) {
    BOOST_CHECK(graph_optimal_result[i] == graph[i]);
  }
  delete operation_graph;
  operation_graph = NULL;

  std::cout << "Test default knn graph kernel with select statements..."
            << std::endl;
  for (std::string &platformName : (*parameters)["PLATFORMS"].keys()) {
    json::Node &platformNode = (*parameters)["PLATFORMS"][platformName];
    for (std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      const std::string &kernelName = "connectNeighbors";
      json::Node &kernelNode = deviceNode["KERNELS"][kernelName];
      kernelNode.replaceIDAttr("USE_SELECT", true);
      kernelNode.replaceIDAttr("VERBOSE", true);
    }
  }
  operation_graph =
      new sgpp::datadriven::DensityOCLMultiPlatform::
        OperationCreateGraphOCLSingleDevice<double>(dataset, 2, manager,
                                                    parameters.get(), k, 0, 0);
  // Test graph kernel
  operation_graph->create_graph(graph);
  for (size_t i = 0; i < dataset.getNrows() * k; ++i) {
    BOOST_CHECK(graph_optimal_result[i] == graph[i]);
  }

  delete operation_graph;
  operation_graph = NULL;

  std::cout << "Test default knn graph kernel with local memory..."
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
      new sgpp::datadriven::DensityOCLMultiPlatform::
        OperationCreateGraphOCLSingleDevice<double>(dataset, 2, manager,
                                                    parameters.get(), k, 0, 0);
  // Test graph kernel
  operation_graph->create_graph(graph);
  for (size_t i = 0; i < dataset.getNrows() * k; ++i) {
    BOOST_CHECK(graph_optimal_result[i] == graph[i]);
  }
  delete operation_graph;
  operation_graph = NULL;

  std::cout << "Test default knn graph kernel with local memory and select statements..."
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
      new sgpp::datadriven::DensityOCLMultiPlatform::
        OperationCreateGraphOCLSingleDevice<double>(dataset, 2, manager,
                                                    parameters.get(), k, 0, 0);
  // Test graph kernel
  operation_graph->create_graph(graph);
  for (size_t i = 0; i < dataset.getNrows() * k; ++i) {
    BOOST_CHECK(graph_optimal_result[i] == graph[i]);
  }
  delete operation_graph;
}

BOOST_AUTO_TEST_CASE(KNNPruneGraphOpenCL)  {
  // Load input
  std::vector<int> graph;
  std::ifstream graph_in( "/media/DatenPartition/Programmieren/ipvs/SGpp/datadriven/tests/datasets/graph_erg_dim2_depth11.txt");
  if (graph_in) {
    double value;
    while (graph_in >> value)
      graph.push_back(value);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("knn graph result file is missing!"));
  }

  // Load optimal results for comparison
  std::vector<int> graph_optimal_result;
  std::ifstream graph_result( "/media/DatenPartition/Programmieren/ipvs/SGpp/datadriven/tests/datasets/graph_pruned_erg_dim2_depth11.txt");
  if (graph_result) {
    double value;
    while (graph_result >> value)
      graph_optimal_result.push_back(value);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("knn pruned graph result file is missing!"));
  }
  std::cerr << "Loaded graphs" << std::endl;

  // Create grid for test scenario
  sgpp::base::Grid *grid = sgpp::base::Grid::createLinearGrid(2);
  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(11);
  size_t gridsize = grid->getStorage().getSize();
    std::cerr << "start " << std::endl;

  // Load optimal results for comparison
  sgpp::base::DataVector alpha(gridsize);
  std::ifstream alpha_result( "/media/DatenPartition/Programmieren/ipvs/SGpp/datadriven/tests/datasets/alpha_erg_dim2_depth11.txt");
  if (alpha_result) {
    int counter = 0;
    std::cerr << "start counter " << counter << std::endl;
    double value;
    while (alpha_result >> value) {
      alpha[counter] = value;
      counter++;
    }
    std::cerr << "End counter " << counter << std::endl;
    BOOST_CHECK(counter == gridsize);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("knn pruned graph result file is missing!"));
  }
  alpha_result.close();
  std::cerr << "Loaded alphas" << std::endl;

  // Create OCL configuration
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();
  sgpp::datadriven::DensityOCLMultiPlatform::
      OperationPruneGraphOCL::load_default_parameters(parameters.get());

  // Create OpenCL Manager
  std::shared_ptr<sgpp::base::OCLManagerMultiPlatform> manager;
  manager = std::make_shared<sgpp::base::OCLManagerMultiPlatform>(true);

  // Load dataset for test scenario
  sgpp::datadriven::Dataset data =
      sgpp::datadriven::ARFFTools::
      readARFF("/media/DatenPartition/Programmieren/ipvs/SGpp/datadriven/tests/datasets/clustering_testdataset_dim2.arff");
  sgpp::base::DataMatrix& dataset = data.getData();

  // Create operation
  sgpp::datadriven::DensityOCLMultiPlatform::OperationPruneGraphOCL* operation_prune =
      new sgpp::datadriven::DensityOCLMultiPlatform::
        OperationPruneGraphOCLMultiPlatform<double>(*grid, alpha, dataset, 2, manager,
                                                    parameters.get(), 0.2, 8, 0, 0);

  std::cout << "Test knn prune kernel ..." << std::endl;
  operation_prune->prune_graph(graph);
  for (size_t i = 0; i < gridsize; ++i) {
    BOOST_CHECK(graph[i] == graph_optimal_result[i]);
  }
  delete operation_prune;
}

BOOST_AUTO_TEST_SUITE_END()
#endif
