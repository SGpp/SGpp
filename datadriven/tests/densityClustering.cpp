// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#if USE_OCL == 1

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <zlib.h>

#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>
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

BOOST_AUTO_TEST_SUITE(TestClusteringOpenCL)

BOOST_AUTO_TEST_CASE(DensityMultiplicationOpenCL) {
  std::vector<std::string> fileNames = {"datadriven/tests/data/clustering_testdataset_dim2.arff"};

  // Create OCL configuration
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();
  sgpp::datadriven::DensityOCLMultiPlatform::
      OperationDensityOCL::load_default_parameters(parameters.get());
  std::vector<std::reference_wrapper<json::Node>> deviceNodes = parameters->getAllDeviceNodes();

  // Create OpenCL Manager
  std::shared_ptr<sgpp::base::OCLManagerMultiPlatform> manager;
  manager = std::make_shared<sgpp::base::OCLManagerMultiPlatform>(true);

  // Create grid for test scenario
  std::unique_ptr<sgpp::base::Grid> grid = sgpp::base::Grid::createLinearGrid(2);
  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(11);
  size_t gridsize = grid->getStorage().getSize();

  // Create vectors for multiplication
  sgpp::base::DataVector alpha(gridsize);
  sgpp::base::DataVector result(gridsize);
  alpha.setAll(1.0);

  // Create operation
  sgpp::datadriven::DensityOCLMultiPlatform::OperationDensityOCL* mult_operation =
      new sgpp::datadriven::DensityOCLMultiPlatform::
      OperationDensityOCLMultiPlatform<double>(*grid, 2, manager, parameters.get(), 0.0001,
                                               0, 0);

  // Execute multiplication
  mult_operation->mult(alpha, result);

  // Load correct results for comparison
  std::vector<double> mult_optimal_result;
  std::ifstream mult_in("datadriven/tests/mult_erg_mult_erg_dim2_depth11.txt");
  if (mult_in) {
    double value;
    while (mult_in >> value)
      mult_optimal_result.push_back(value);

    // Compare results with optimal results
    for (size_t i = 0; i < gridsize; ++i) {
      BOOST_CHECK(mult_optimal_result[i] == result[i]);
    }
  } else {
    std::cerr << "Density multiplication result file is missing!" << std::endl;
    std::cerr << "Cannot execute this unit test without these results..." << std::endl;
  }
}
BOOST_AUTO_TEST_SUITE_END()

#endif
