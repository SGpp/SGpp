// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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
#include <sgpp/datadriven/algorithm/GridFactory.hpp>
#include <sgpp/globaldef.hpp>

#include <string>
#include <set>
#include <vector>

BOOST_AUTO_TEST_SUITE(dBMatOffline_test)

BOOST_AUTO_TEST_CASE(testReadWriteOrthoAdapt) {
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 2;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  sgpp::base::AdaptivityConfiguration adaptivityConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 0.1;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;

  sgpp::datadriven::GridFactory gridFactory;
  std::unique_ptr<sgpp::base::Grid> grid = std::unique_ptr<sgpp::base::Grid>{
      gridFactory.createGrid(gridConfig, std::set<std::set<size_t>>())};

  auto offline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(
          gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig)};
  offline->buildMatrix(grid.get(), regularizationConfig);
  offline->decomposeMatrix(regularizationConfig, densityEstimationConfig);
  std::string filename = "test.dbmat";
  offline->store(filename);
  auto newOffline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildFromFile(filename)};
  std::remove(filename.c_str());

  /**
   * Check matrices
   */
  // lhsMatrix of parent object
  auto& oldMatrix = offline->getDecomposedMatrix();

  std::cout << "Got old matrix decom" << std::endl;
  auto& newMatrix = newOffline->getDecomposedMatrix();

  std::cout << "Got decompositions" << std::endl;

  BOOST_CHECK_EQUAL(oldMatrix.getSize(), newMatrix.getSize());

  for (size_t i = 0; i < newMatrix.getSize(); i++) {
    BOOST_CHECK_CLOSE(newMatrix[i], oldMatrix[i], 1e-4);
  }

  std::cout << "Are close" << std::endl;

  // q_ortho_matrix_
  sgpp::datadriven::DBMatOfflineOrthoAdapt* child =
      static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt*>(&*offline);
  auto& oldMatrixQ = child->getQ();
  sgpp::datadriven::DBMatOfflineOrthoAdapt* newChild =
      static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt*>(&*newOffline);
  auto& newMatrixQ = newChild->getQ();

  std::cout << "Got Q" << std::endl;

  BOOST_CHECK_EQUAL(oldMatrixQ.getSize(), newMatrixQ.getSize());

  for (size_t i = 0; i < newMatrixQ.getSize(); i++) {
    BOOST_CHECK_CLOSE(newMatrixQ[i], oldMatrixQ[i], 1e-4);
  }

  // t_tridiag_inv
  auto& oldMatrixTinv = child->getTinv();
  auto& newMatrixTinv = newChild->getTinv();

  std::cout << "Got T" << std::endl;

  BOOST_CHECK_EQUAL(oldMatrixTinv.getSize(), newMatrixTinv.getSize());

  for (size_t i = 0; i < newMatrixTinv.getSize(); i++) {
    BOOST_CHECK_CLOSE(newMatrixTinv[i], oldMatrixTinv[i], 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(testReadWriteCholesky) {
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 2;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  sgpp::base::AdaptivityConfiguration adaptivityConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 0.1;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::Chol;

  sgpp::datadriven::GridFactory gridFactory;
  std::unique_ptr<sgpp::base::Grid> grid = std::unique_ptr<sgpp::base::Grid>{
      gridFactory.createGrid(gridConfig, std::set<std::set<size_t>>())};

  auto offline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(
          gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig)};
  offline->buildMatrix(grid.get(), regularizationConfig);
  offline->decomposeMatrix(regularizationConfig, densityEstimationConfig);

  std::string filename = "test.dbmat";
  offline->store(filename);
  auto newOffline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildFromFile(filename)};
  std::remove(filename.c_str());

  /**
   * Check matrices
   */
  auto& oldMatrix = offline->getDecomposedMatrix();
  auto& newMatrix = newOffline->getDecomposedMatrix();

  BOOST_CHECK_EQUAL(oldMatrix.getSize(), newMatrix.getSize());

  for (size_t i = 0; i < newMatrix.getSize(); i++) {
    BOOST_CHECK_CLOSE(newMatrix[i], oldMatrix[i], 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(testReadWriteEigen) {
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 2;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  sgpp::base::AdaptivityConfiguration adaptivityConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 0.1;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::Eigen;

  sgpp::datadriven::GridFactory gridFactory;
  std::unique_ptr<sgpp::base::Grid> grid = std::unique_ptr<sgpp::base::Grid>{
      gridFactory.createGrid(gridConfig, std::set<std::set<size_t>>())};

  auto offline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(
          gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig)};
  offline->buildMatrix(grid.get(), regularizationConfig);
  offline->decomposeMatrix(regularizationConfig, densityEstimationConfig);

  std::string filename = "test.dbmat";
  offline->store(filename);
  auto newOffline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildFromFile(filename)};
  std::remove(filename.c_str());

  /**
   * Check matrices
   */
  auto& oldMatrix = offline->getDecomposedMatrix();
  auto& newMatrix = newOffline->getDecomposedMatrix();

  BOOST_CHECK_EQUAL(oldMatrix.getSize(), newMatrix.getSize());

  for (size_t i = 0; i < newMatrix.getSize(); i++) {
    BOOST_CHECK_CLOSE(newMatrix[i], oldMatrix[i], 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(testReadWriteLU) {
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 2;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  sgpp::base::AdaptivityConfiguration adaptivityConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 0.1;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::LU;

  sgpp::datadriven::GridFactory gridFactory;
  std::unique_ptr<sgpp::base::Grid> grid = std::unique_ptr<sgpp::base::Grid>{
      gridFactory.createGrid(gridConfig, std::set<std::set<size_t>>())};

  auto offline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(
          gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig)};
  offline->buildMatrix(grid.get(), regularizationConfig);
  offline->decomposeMatrix(regularizationConfig, densityEstimationConfig);

  std::string filename = "test.dbmat";
  offline->store(filename);
  auto newOffline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
      sgpp::datadriven::DBMatOfflineFactory::buildFromFile(filename)};
  std::remove(filename.c_str());

  /**
   * Check matrices
   */
  auto& oldMatrix = offline->getDecomposedMatrix();
  auto& newMatrix = newOffline->getDecomposedMatrix();

  BOOST_CHECK_EQUAL(oldMatrix.getSize(), newMatrix.getSize());

  for (size_t i = 0; i < newMatrix.getSize(); i++) {
    BOOST_CHECK_CLOSE(newMatrix[i], oldMatrix[i], 1e-4);
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* USE_GSL */
