/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * test_DBMatOffline.cpp
 *
 *  Created on: Apr 8, 2017
 *      Author: Michael Lettrich
 */
#ifdef USE_GSL
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineLU.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/globaldef.hpp>

#include <string>

BOOST_AUTO_TEST_SUITE(dBMatOffline_test)

BOOST_AUTO_TEST_CASE(testReadWriteCholesky) {
  sgpp::datadriven::DBMatDensityConfiguration config;
  config.grid_dim_ = 2;
  config.grid_level_ = 3;
  config.grid_type_ = sgpp::base::GridType::Linear;
  config.regularization_ = sgpp::datadriven::RegularizationType::Identity;
  config.lambda_ = 0.1;
  config.decomp_type_ = sgpp::datadriven::DBMatDecompostionType::Chol;

  auto offline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(config)};
  offline->buildMatrix();
  offline->decomposeMatrix();

  std::string filename = "test.dbmat";
  offline->store(filename);
  auto newOffline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildFromFile(filename)};

  auto newConfig = newOffline->getConfig();

  /**
   * Check Configuration
   */
  BOOST_CHECK_EQUAL(config.grid_dim_, newConfig.grid_dim_);
  BOOST_CHECK_EQUAL(config.grid_level_, newConfig.grid_level_);
  BOOST_CHECK_EQUAL(static_cast<int>(config.grid_type_), static_cast<int>(newConfig.grid_type_));
  BOOST_CHECK_EQUAL(static_cast<int>(config.regularization_),
                    static_cast<int>(newConfig.regularization_));
  BOOST_CHECK_CLOSE(config.lambda_, newConfig.lambda_, 10e-5);
  BOOST_CHECK_EQUAL(static_cast<int>(config.decomp_type_),
                    static_cast<int>(newConfig.decomp_type_));

  /**
   * Check matrices
   */
  auto& oldMatrix = offline->getDecomposedMatrix();
  auto& newMatrix = newOffline->getDecomposedMatrix();

  BOOST_CHECK_EQUAL(oldMatrix.getSize(), newMatrix.getSize());

  for (size_t i = 0; i < newMatrix.getSize(); i++) {
    BOOST_CHECK_CLOSE(newMatrix[i], oldMatrix[i], 10e-5);
  }
}

BOOST_AUTO_TEST_CASE(testReadWriteEigen) {
  sgpp::datadriven::DBMatDensityConfiguration config;
  config.grid_dim_ = 2;
  config.grid_level_ = 3;
  config.grid_type_ = sgpp::base::GridType::Linear;
  config.regularization_ = sgpp::datadriven::RegularizationType::Identity;
  config.lambda_ = 0.1;
  config.decomp_type_ = sgpp::datadriven::DBMatDecompostionType::Eigen;

  auto offline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(config)};
  offline->buildMatrix();
  offline->decomposeMatrix();

  std::string filename = "test.dbmat";
  offline->store(filename);
  auto newOffline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildFromFile(filename)};

  auto newConfig = newOffline->getConfig();

  /**
   * Check Configuration
   */
  BOOST_CHECK_EQUAL(config.grid_dim_, newConfig.grid_dim_);
  BOOST_CHECK_EQUAL(config.grid_level_, newConfig.grid_level_);
  BOOST_CHECK_EQUAL(static_cast<int>(config.grid_type_), static_cast<int>(newConfig.grid_type_));
  BOOST_CHECK_EQUAL(static_cast<int>(config.regularization_),
                    static_cast<int>(newConfig.regularization_));
  BOOST_CHECK_CLOSE(config.lambda_, newConfig.lambda_, 10e-5);
  BOOST_CHECK_EQUAL(static_cast<int>(config.decomp_type_),
                    static_cast<int>(newConfig.decomp_type_));

  /**
   * Check matrices
   */
  auto& oldMatrix = offline->getDecomposedMatrix();
  auto& newMatrix = newOffline->getDecomposedMatrix();

  BOOST_CHECK_EQUAL(oldMatrix.getSize(), newMatrix.getSize());

  for (size_t i = 0; i < newMatrix.getSize(); i++) {
    BOOST_CHECK_CLOSE(newMatrix[i], oldMatrix[i], 10e-5);
  }
}

BOOST_AUTO_TEST_CASE(testReadWriteLU) {
  sgpp::datadriven::DBMatDensityConfiguration config;
  config.grid_dim_ = 2;
  config.grid_level_ = 3;
  config.grid_type_ = sgpp::base::GridType::Linear;
  config.regularization_ = sgpp::datadriven::RegularizationType::Identity;
  config.lambda_ = 0.1;
  config.decomp_type_ = sgpp::datadriven::DBMatDecompostionType::LU;

  auto offline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(config)};
  offline->buildMatrix();
  offline->decomposeMatrix();

  std::string filename = "test.dbmat";
  offline->store(filename);
  auto newOffline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildFromFile(filename)};

  auto newConfig = newOffline->getConfig();

  /**
   * Check Configuration
   */
  BOOST_CHECK_EQUAL(config.grid_dim_, newConfig.grid_dim_);
  BOOST_CHECK_EQUAL(config.grid_level_, newConfig.grid_level_);
  BOOST_CHECK_EQUAL(static_cast<int>(config.grid_type_), static_cast<int>(newConfig.grid_type_));
  BOOST_CHECK_EQUAL(static_cast<int>(config.regularization_),
                    static_cast<int>(newConfig.regularization_));
  BOOST_CHECK_CLOSE(config.lambda_, newConfig.lambda_, 10e-5);
  BOOST_CHECK_EQUAL(static_cast<int>(config.decomp_type_),
                    static_cast<int>(newConfig.decomp_type_));

  /**
   * Check matrices
   */
  auto& oldMatrix = offline->getDecomposedMatrix();
  auto& newMatrix = newOffline->getDecomposedMatrix();

  BOOST_CHECK_EQUAL(oldMatrix.getSize(), newMatrix.getSize());

  for (size_t i = 0; i < newMatrix.getSize(); i++) {
    BOOST_CHECK_CLOSE(newMatrix[i], oldMatrix[i], 10e-5);
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif /*USE_GSL*/
