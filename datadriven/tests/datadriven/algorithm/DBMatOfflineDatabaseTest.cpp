// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using sgpp::base::AdaptivityConfiguration;
using sgpp::base::GeneralGridConfiguration;
using sgpp::datadriven::DBMatDatabase;
using sgpp::datadriven::DensityEstimationConfiguration;
using sgpp::datadriven::RegularizationConfiguration;

BOOST_AUTO_TEST_SUITE(DBMatDatabaseTest)

std::string createEmptyDatabase() {
  // Create an empty database
  // The file will be placed in the cwd temporarily to ensure cross plattform compatibility
  std::string filepath = "tmpdatabase";
  std::ofstream stream(filepath);
  stream << "{" << std::endl << "\"database\" : []" << std::endl << "}" << std::endl;
  stream.close();
  return filepath;
}

void removeDatabase(std::string& path) { remove(path.c_str()); }

void initializeStandardConfiguration(GeneralGridConfiguration& gridConfig,
                                     AdaptivityConfiguration& adaptConfig,
                                     RegularizationConfiguration& regularizationConfig,
                                     DensityEstimationConfiguration& densityEstimationConfig) {
  gridConfig.dim_ = 5;
  gridConfig.level_ = 3;
  regularizationConfig.lambda_ = 1e-2;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;
  densityEstimationConfig.type_ = sgpp::datadriven::DensityEstimationType::Decomposition;
}

BOOST_AUTO_TEST_CASE(TestInitialization) {
  // Create a empty database
  std::string filepath = createEmptyDatabase();
  DBMatDatabase database(filepath);
  removeDatabase(filepath);
}

BOOST_AUTO_TEST_CASE(TestPutgMatrix) {
  // Test to put a certain datamatrix
  std::string filepath = createEmptyDatabase();
  DBMatDatabase database(filepath);
  sgpp::base::RegularGridConfiguration gridConfig;
  AdaptivityConfiguration adaptConfig;
  RegularizationConfiguration regularizationConfig;
  DensityEstimationConfiguration densityEstimationConfig;
  initializeStandardConfiguration(gridConfig, adaptConfig, regularizationConfig,
                                  densityEstimationConfig);
  database.putDataMatrix(gridConfig, adaptConfig, regularizationConfig,
                         densityEstimationConfig, "testfilepath");

  // Assert that the database holds a matrix
  BOOST_CHECK(database.hasDataMatrix(gridConfig, adaptConfig, regularizationConfig,
                                     densityEstimationConfig));
  removeDatabase(filepath);
}

BOOST_AUTO_TEST_CASE(TestGetMatrix) {
  // Test to put a certain datamatrix
  std::string filepath = createEmptyDatabase();
  DBMatDatabase database(filepath);
  sgpp::base::RegularGridConfiguration gridConfig;
  AdaptivityConfiguration adaptConfig;
  RegularizationConfiguration regularizationConfig;
  DensityEstimationConfiguration densityEstimationConfig;

  std::string path = "testfilepathGetMatrix";

  initializeStandardConfiguration(gridConfig, adaptConfig, regularizationConfig,
                                  densityEstimationConfig);
  database.putDataMatrix(gridConfig, adaptConfig, regularizationConfig,
                         densityEstimationConfig, path);

  // Assert that the database holds a matrix
  BOOST_CHECK(database.hasDataMatrix(gridConfig, adaptConfig, regularizationConfig,
                                     densityEstimationConfig));

  std::string& pathFromDatabase = database.getDataMatrix(
      gridConfig, adaptConfig, regularizationConfig, densityEstimationConfig);

  // Assert the retrieved string equals the string passed
  BOOST_CHECK(path.compare(pathFromDatabase) == 0);

  // Modify the configuration and assert the database does not hold the string
  regularizationConfig.lambda_ *= 2.0;
  BOOST_CHECK(!database.hasDataMatrix(gridConfig, adaptConfig, regularizationConfig,
                                      densityEstimationConfig));
  removeDatabase(filepath);
}

BOOST_AUTO_TEST_CASE(TestOverwriteMatrix) {
  // Test to put a certain datamatrix
  std::string filepath = createEmptyDatabase();
  DBMatDatabase database(filepath);
  sgpp::base::RegularGridConfiguration gridConfig;
  AdaptivityConfiguration adaptConfig;
  RegularizationConfiguration regularizationConfig;
  DensityEstimationConfiguration densityEstimationConfig;

  std::string path = "testfilepathGetMatrix";

  initializeStandardConfiguration(gridConfig, adaptConfig, regularizationConfig,
                                  densityEstimationConfig);
  database.putDataMatrix(gridConfig, adaptConfig, regularizationConfig,
                         densityEstimationConfig, path);

  // Assert that the database holds a matrix
  BOOST_CHECK(database.hasDataMatrix(gridConfig, adaptConfig, regularizationConfig,
                                     densityEstimationConfig));

  std::string& pathFromDatabase = database.getDataMatrix(
      gridConfig, adaptConfig, regularizationConfig, densityEstimationConfig);

  // Assert the retrieved string equals the string passed
  BOOST_CHECK(path.compare(pathFromDatabase) == 0);

  // Change the path associated with the configuration
  std::string newpath = "newpathGetMatrix";
  database.putDataMatrix(gridConfig, adaptConfig, regularizationConfig,
                         densityEstimationConfig, newpath, true);

  pathFromDatabase = database.getDataMatrix(gridConfig, adaptConfig, regularizationConfig,
                                            densityEstimationConfig);

  // Assert the path was overwritten
  BOOST_CHECK(newpath.compare(pathFromDatabase) == 0);
  removeDatabase(filepath);
}

BOOST_AUTO_TEST_CASE(TestBaseMatrix) {
  std::string filepath = createEmptyDatabase();
  DBMatDatabase database(filepath);

  sgpp::base::AdaptivityConfiguration adaptConfig;
  adaptConfig.numRefinements_ = 0;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.lambda_ = 0;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.type_ = sgpp::datadriven::DensityEstimationType::Decomposition;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;

  sgpp::base::GeneralGridConfiguration baseGridConfig1;
  baseGridConfig1.generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
  baseGridConfig1.type_ = sgpp::base::GridType::Linear;
  baseGridConfig1.levelVector_ = std::vector<size_t>{2, 3, 3};
  baseGridConfig1.dim_ = 3;

  std::string testPath1 = "testFilePathBaseMatrix1";

  database.putDataMatrix(baseGridConfig1, adaptConfig, regularizationConfig,
                         densityEstimationConfig, testPath1);

  sgpp::base::GeneralGridConfiguration baseGridConfig2;
  baseGridConfig2.generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
  baseGridConfig2.type_ = sgpp::base::GridType::Linear;
  baseGridConfig2.levelVector_ = std::vector<size_t>{2, 2, 3};
  baseGridConfig2.dim_ = 3;

  std::string testPath2 = "testFilePathBaseMatrix2";

  database.putDataMatrix(baseGridConfig2, adaptConfig, regularizationConfig,
                         densityEstimationConfig, testPath2);

  // Test case: Permutation
  sgpp::base::GeneralGridConfiguration testGridConfig1;
  testGridConfig1.generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
  testGridConfig1.type_ = sgpp::base::GridType::Linear;
  testGridConfig1.levelVector_ = std::vector<size_t>{3, 2, 3};
  testGridConfig1.dim_ = 3;

  BOOST_CHECK(database.hasBaseDataMatrix(testGridConfig1, adaptConfig, regularizationConfig,
                                         densityEstimationConfig));

  sgpp::base::GeneralGridConfiguration obtainedBaseGridConfig;

  std::string obtainedPath =
      database.getBaseDataMatrix(testGridConfig1, adaptConfig, regularizationConfig,
                                 densityEstimationConfig, obtainedBaseGridConfig);

  BOOST_CHECK(obtainedPath.compare(testPath1) == 0);
  BOOST_CHECK(baseGridConfig1.levelVector_ == obtainedBaseGridConfig.levelVector_ &&
              baseGridConfig1.dim_ == obtainedBaseGridConfig.dim_);

  // Test case: Blow-up
  sgpp::base::GeneralGridConfiguration testGridConfig2;
  testGridConfig2.generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
  testGridConfig2.type_ = sgpp::base::GridType::Linear;
  testGridConfig2.levelVector_ = std::vector<size_t>{2, 2, 3, 1, 1, 1, 1, 1};
  testGridConfig2.dim_ = 8;

  BOOST_CHECK(database.hasBaseDataMatrix(testGridConfig2, adaptConfig, regularizationConfig,
                                         densityEstimationConfig));

  obtainedPath = database.getBaseDataMatrix(testGridConfig2, adaptConfig, regularizationConfig,
                                            densityEstimationConfig, obtainedBaseGridConfig);

  BOOST_CHECK(obtainedPath.compare(testPath2) == 0);
  BOOST_CHECK(baseGridConfig2.levelVector_ == obtainedBaseGridConfig.levelVector_ &&
              baseGridConfig2.dim_ == obtainedBaseGridConfig.dim_);

  // Test case: Permutation and blow-up
  sgpp::base::GeneralGridConfiguration testGridConfig3;
  testGridConfig3.generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
  testGridConfig3.type_ = sgpp::base::GridType::Linear;
  testGridConfig3.levelVector_ = std::vector<size_t>{3, 3, 2, 1, 1, 1, 1, 1};
  testGridConfig3.dim_ = 8;

  BOOST_CHECK(database.hasBaseDataMatrix(testGridConfig3, adaptConfig, regularizationConfig,
                                         densityEstimationConfig));

  obtainedPath = database.getBaseDataMatrix(testGridConfig3, adaptConfig, regularizationConfig,
                                            densityEstimationConfig, obtainedBaseGridConfig);

  BOOST_CHECK(obtainedPath.compare(testPath1) == 0);
  BOOST_CHECK(baseGridConfig1.levelVector_ == obtainedBaseGridConfig.levelVector_ &&
              baseGridConfig1.dim_ == obtainedBaseGridConfig.dim_);

  // Test case: No base
  sgpp::base::GeneralGridConfiguration testGridConfig4;
  testGridConfig3.generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
  testGridConfig3.type_ = sgpp::base::GridType::Linear;
  testGridConfig3.levelVector_ = std::vector<size_t>{3, 3, 3, 1, 1, 1, 1, 1};
  testGridConfig3.dim_ = 8;

  BOOST_CHECK(!database.hasBaseDataMatrix(testGridConfig3, adaptConfig, regularizationConfig,
                                          densityEstimationConfig));
  removeDatabase(filepath);
}

BOOST_AUTO_TEST_SUITE_END()
