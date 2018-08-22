// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <iostream>
#include <cstdio>
#include <string>

using sgpp::base::GeneralGridConfiguration;
using sgpp::base::AdaptivityConfiguration;
using sgpp::datadriven::RegularizationConfiguration;
using sgpp::datadriven::DensityEstimationConfiguration;
using sgpp::datadriven::DBMatDatabase;

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

void removeDatabase(std::string& path) {
  remove(path.c_str());
}

void initializeStandardConfiguration(
    GeneralGridConfiguration &gridConfig,
    AdaptivityConfiguration& adaptivityConfig,
    RegularizationConfiguration& regularizationConfig,
    DensityEstimationConfiguration& densityEstimationConfig
    ) {
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
  AdaptivityConfiguration adaptivityConfig;
  RegularizationConfiguration regularizationConfig;
  DensityEstimationConfiguration densityEstimationConfig;
  initializeStandardConfiguration(gridConfig, adaptivityConfig, regularizationConfig,
      densityEstimationConfig);
  database.putDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
      densityEstimationConfig, "testfilepath");

  // Assert that the database holds a matrix
  BOOST_CHECK(database.hasDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
      densityEstimationConfig));
  removeDatabase(filepath);
}

BOOST_AUTO_TEST_CASE(TestGetMatrix) {
  // Test to put a certain datamatrix
  std::string filepath = createEmptyDatabase();
  DBMatDatabase database(filepath);
  sgpp::base::RegularGridConfiguration gridConfig;
  AdaptivityConfiguration adaptivityConfig;
  RegularizationConfiguration regularizationConfig;
  DensityEstimationConfiguration densityEstimationConfig;

  std::string path = "testfilepathGetMatrix";

  initializeStandardConfiguration(gridConfig, adaptivityConfig, regularizationConfig,
      densityEstimationConfig);
  database.putDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
      densityEstimationConfig, path);

  // Assert that the database holds a matrix
  BOOST_CHECK(database.hasDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
      densityEstimationConfig));

  std::string& pathFromDatabase = database.getDataMatrix(gridConfig, adaptivityConfig,
      regularizationConfig, densityEstimationConfig);

  // Assert the retrieved string equals the string passed
  BOOST_CHECK(path.compare(pathFromDatabase) == 0);

  // Modify the configuration and assert the database does not hold the string
  regularizationConfig.lambda_ *= 2.0;
  BOOST_CHECK(!database.hasDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
      densityEstimationConfig));
  removeDatabase(filepath);
}

BOOST_AUTO_TEST_CASE(TestOverwriteMatrix) {
  // Test to put a certain datamatrix
  std::string filepath = createEmptyDatabase();
  DBMatDatabase database(filepath);
  sgpp::base::RegularGridConfiguration gridConfig;
  AdaptivityConfiguration adaptivityConfig;
  RegularizationConfiguration regularizationConfig;
  DensityEstimationConfiguration densityEstimationConfig;

  std::string path = "testfilepathGetMatrix";

  initializeStandardConfiguration(gridConfig, adaptivityConfig, regularizationConfig,
      densityEstimationConfig);
  database.putDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
      densityEstimationConfig, path);

  // Assert that the database holds a matrix
  BOOST_CHECK(database.hasDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
      densityEstimationConfig));

  std::string& pathFromDatabase = database.getDataMatrix(gridConfig, adaptivityConfig,
      regularizationConfig, densityEstimationConfig);

  // Assert the retrieved string equals the string passed
  BOOST_CHECK(path.compare(pathFromDatabase) == 0);

  // Change the path associated with the configuration
  std::string newpath = "newpathGetMatrix";
  database.putDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
       densityEstimationConfig, newpath, true);

  pathFromDatabase = database.getDataMatrix(gridConfig, adaptivityConfig,
        regularizationConfig, densityEstimationConfig);

  // Assert the path was overwritten
  BOOST_CHECK(newpath.compare(pathFromDatabase) == 0);
  removeDatabase(filepath);
}

BOOST_AUTO_TEST_SUITE_END()
