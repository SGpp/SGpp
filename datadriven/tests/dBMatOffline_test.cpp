/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * test_DBMatOffline.cpp
 *
 * Created on: Apr 8, 2017
 *      Author: Michael Lettrich
 */
#ifdef USE_GSL

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineLU.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/globaldef.hpp>

#include <string>

BOOST_AUTO_TEST_SUITE(dBMatOffline_test)

BOOST_AUTO_TEST_CASE(testReadWriteOrthoAdapt) {
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 2;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  sgpp::base::AdpativityConfiguration adaptivityConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 0.1;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;

  auto offline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(gridConfig,
                                                                adaptivityConfig,
                                                                regularizationConfig,
                                                                densityEstimationConfig)};
  offline->buildMatrix();
  offline->decomposeMatrix();
  std::string filename = "test.dbmat";
  offline->store(filename);
  auto newOffline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildFromFile(filename)};
  std::remove(filename.c_str());
  auto newGridConfig = newOffline->getGridConfig();
  auto newRegularizationConfig = newOffline->getRegularizationConfig();
  auto newDensityEstimationConfig = newOffline->getDensityEstimationConfig();

  /**
   * Check Configuration
   */
  BOOST_CHECK_EQUAL(gridConfig.dim_, newGridConfig.dim_);
  BOOST_CHECK_EQUAL(gridConfig.level_, newGridConfig.level_);
  BOOST_CHECK_EQUAL(static_cast<int>(gridConfig.type_), static_cast<int>(newGridConfig.type_));
  BOOST_CHECK_EQUAL(static_cast<int>(regularizationConfig.type_),
                    static_cast<int>(newRegularizationConfig.type_));
  BOOST_CHECK_CLOSE(regularizationConfig.lambda_, newRegularizationConfig.lambda_, 10e-5);
  BOOST_CHECK_EQUAL(static_cast<int>(densityEstimationConfig.decomposition_),
                    static_cast<int>(newDensityEstimationConfig.decomposition_));
  /**
   * Check matrices
   */
  // lhsMatrix of parent object
  auto& oldMatrix = offline->getDecomposedMatrix();
  auto& newMatrix = newOffline->getDecomposedMatrix();

  BOOST_CHECK_EQUAL(oldMatrix.getSize(), newMatrix.getSize());

  for (size_t i = 0; i < newMatrix.getSize(); i++) {
    BOOST_CHECK_CLOSE(newMatrix[i], oldMatrix[i], 10e-5);
  }

  // q_ortho_matrix_
  sgpp::datadriven::DBMatOfflineOrthoAdapt* child =
      static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt*>(&*offline);
  auto& oldMatrixQ = child->getQ();
  sgpp::datadriven::DBMatOfflineOrthoAdapt* newChild =
      static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt*>(&*newOffline);
  auto& newMatrixQ = newChild->getQ();

  BOOST_CHECK_EQUAL(oldMatrixQ.getSize(), newMatrixQ.getSize());

  for (size_t i = 0; i < newMatrixQ.getSize(); i++) {
    BOOST_CHECK_CLOSE(newMatrixQ[i], oldMatrixQ[i], 10e-5);
  }

  // t_tridiag_inv
  auto& oldMatrixTinv = child->getTinv();
  auto& newMatrixTinv = newChild->getTinv();

  BOOST_CHECK_EQUAL(oldMatrixTinv.getSize(), newMatrixTinv.getSize());

  for (size_t i = 0; i < newMatrixTinv.getSize(); i++) {
    BOOST_CHECK_CLOSE(newMatrixTinv[i], oldMatrixTinv[i], 10e-5);
  }
}

BOOST_AUTO_TEST_CASE(testReadWriteCholesky) {
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 2;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  sgpp::base::AdpativityConfiguration adaptivityConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 0.1;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::Chol;

  auto offline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(gridConfig,
                                                                adaptivityConfig,
                                                                regularizationConfig,
                                                                densityEstimationConfig)};
  offline->buildMatrix();
  offline->decomposeMatrix();

  std::string filename = "test.dbmat";
  offline->store(filename);
  auto newOffline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildFromFile(filename)};
  std::remove(filename.c_str());
  auto newGridConfig = newOffline->getGridConfig();
  auto newRegularizationConfig = newOffline->getRegularizationConfig();
  auto newDensityEstimationConfig = newOffline->getDensityEstimationConfig();
  /**
   * Check Configuration
   */
  BOOST_CHECK_EQUAL(gridConfig.dim_, newGridConfig.dim_);
  BOOST_CHECK_EQUAL(gridConfig.level_, newGridConfig.level_);
  BOOST_CHECK_EQUAL(static_cast<int>(gridConfig.type_), static_cast<int>(newGridConfig.type_));
  BOOST_CHECK_EQUAL(static_cast<int>(regularizationConfig.type_),
                    static_cast<int>(newRegularizationConfig.type_));
  BOOST_CHECK_CLOSE(regularizationConfig.lambda_, newRegularizationConfig.lambda_, 10e-5);
  BOOST_CHECK_EQUAL(static_cast<int>(densityEstimationConfig.decomposition_),
                    static_cast<int>(newDensityEstimationConfig.decomposition_));

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
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 2;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  sgpp::base::AdpativityConfiguration adaptivityConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 0.1;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::Eigen;

  auto offline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(gridConfig,
                                                                adaptivityConfig,
                                                                regularizationConfig,
                                                                densityEstimationConfig)};
  offline->buildMatrix();
  offline->decomposeMatrix();

  std::string filename = "test.dbmat";
  offline->store(filename);
  auto newOffline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildFromFile(filename)};
  std::remove(filename.c_str());
  auto newGridConfig = newOffline->getGridConfig();
  auto newRegularizationConfig = newOffline->getRegularizationConfig();
  auto newDensityEstimationConfig = newOffline->getDensityEstimationConfig();
  /**
   * Check Configuration
   */
  BOOST_CHECK_EQUAL(gridConfig.dim_, newGridConfig.dim_);
  BOOST_CHECK_EQUAL(gridConfig.level_, newGridConfig.level_);
  BOOST_CHECK_EQUAL(static_cast<int>(gridConfig.type_), static_cast<int>(newGridConfig.type_));
  BOOST_CHECK_EQUAL(static_cast<int>(regularizationConfig.type_),
                    static_cast<int>(newRegularizationConfig.type_));
  BOOST_CHECK_CLOSE(regularizationConfig.lambda_, newRegularizationConfig.lambda_, 10e-5);
  BOOST_CHECK_EQUAL(static_cast<int>(densityEstimationConfig.decomposition_),
                    static_cast<int>(newDensityEstimationConfig.decomposition_));

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
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 2;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  sgpp::base::AdpativityConfiguration adaptivityConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 0.1;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::LU;

  auto offline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(gridConfig,
                                                                adaptivityConfig,
                                                                regularizationConfig,
                                                                densityEstimationConfig)};
  offline->buildMatrix();
  offline->decomposeMatrix();

  std::string filename = "test.dbmat";
  offline->store(filename);
  auto newOffline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildFromFile(filename)};
  std::remove(filename.c_str());
  auto newGridConfig = newOffline->getGridConfig();
  auto newRegularizationConfig = newOffline->getRegularizationConfig();
  auto newDensityEstimationConfig = newOffline->getDensityEstimationConfig();
  /**
   * Check Configuration
   */
  BOOST_CHECK_EQUAL(gridConfig.dim_, newGridConfig.dim_);
  BOOST_CHECK_EQUAL(gridConfig.level_, newGridConfig.level_);
  BOOST_CHECK_EQUAL(static_cast<int>(gridConfig.type_), static_cast<int>(newGridConfig.type_));
  BOOST_CHECK_EQUAL(static_cast<int>(regularizationConfig.type_),
                    static_cast<int>(newRegularizationConfig.type_));
  BOOST_CHECK_CLOSE(regularizationConfig.lambda_, newRegularizationConfig.lambda_, 10e-5);
  BOOST_CHECK_EQUAL(static_cast<int>(densityEstimationConfig.decomposition_),
                    static_cast<int>(newDensityEstimationConfig.decomposition_));

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

#endif /* USE_GSL */
