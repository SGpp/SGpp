// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflinePermutable.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/datadriven/configuration/ParallelConfiguration.hpp>

#include <set>
#include <vector>

BOOST_AUTO_TEST_SUITE(DBMatOfflinePermutableTest)

class DBMatOfflineTest : public sgpp::datadriven::DBMatOfflinePermutable {
 public:
  DBMatOfflineTest() {}

  void permuteDecomposition(
      const sgpp::base::GeneralGridConfiguration& baseGridConfig,
      const sgpp::base::GeneralGridConfiguration& desiredGridConfig) override {}

  void decomposeMatrix(
      const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) override {}

  DBMatOffline* clone() const override { return new DBMatOfflineTest(); }

  bool isRefineable() override { return false; }

  const DataMatrix& getUnmodifiedR() override {
    throw sgpp::base::not_implemented_exception(
        "DBMatOfflineTest::getUnmodifiedR() not implemented!");
  }

  const sgpp::datadriven::DataMatrixDistributed& getUnmodifiedRDistributed(
      std::shared_ptr<sgpp::datadriven::BlacsProcessGrid> processGrid,
      const sgpp::datadriven::ParallelConfiguration& parallelConfig) override {
    throw sgpp::base::not_implemented_exception(
        "DBMatOfflineTest::getUnmodifiedRDistributed() not implemented!");
  }

  void updateRegularization(double lambda) override {
    throw sgpp::base::not_implemented_exception(
        "DBMatOfflineTest::updateRegularization() not implemented!");
  }

  void updateRegularizationParallel(
      double lambda, std::shared_ptr<sgpp::datadriven::BlacsProcessGrid> processGrid,
      const sgpp::datadriven::ParallelConfiguration& parallelConfig) override {
    throw sgpp::base::not_implemented_exception(
        "DBMatOfflineTest::updateRegularizationParallel() is not implemented!");
  }

  sgpp::datadriven::MatrixDecompositionType getDecompositionType() override {
    return sgpp::datadriven::MatrixDecompositionType::Chol;
  }

  size_t test_getMatrixIndexForPoint(std::vector<size_t> level, std::vector<size_t> index,
                                     std::vector<size_t> gridLevel,
                                     const std::vector<size_t>& preComputations) {
    return this->getMatrixIndexForPoint(level, index, gridLevel, preComputations);
  }

  std::vector<size_t> test_preComputeMatrixIndexForPoint(std::vector<size_t> level) {
    return this->preComputeMatrixIndexForPoint(level);
  }
};

BOOST_AUTO_TEST_CASE(GetMatrixIndexForPointTest) {
  DBMatOfflineTest testObject;
  std::vector<size_t> pre =
      testObject.test_preComputeMatrixIndexForPoint(std::vector<size_t>{3, 2});
  BOOST_CHECK(testObject.test_getMatrixIndexForPoint(std::vector<size_t>{1, 1},
                                                     std::vector<size_t>{1, 1},
                                                     std::vector<size_t>{3, 2}, pre) == 1);
  BOOST_CHECK(testObject.test_getMatrixIndexForPoint(std::vector<size_t>{2, 1},
                                                     std::vector<size_t>{1, 1},
                                                     std::vector<size_t>{3, 2}, pre) == 2);
  BOOST_CHECK(testObject.test_getMatrixIndexForPoint(std::vector<size_t>{2, 2},
                                                     std::vector<size_t>{1, 3},
                                                     std::vector<size_t>{3, 2}, pre) == 11);
  BOOST_CHECK(testObject.test_getMatrixIndexForPoint(std::vector<size_t>{3, 2},
                                                     std::vector<size_t>{5, 1},
                                                     std::vector<size_t>{3, 2}, pre) == 18);
  BOOST_CHECK(testObject.test_getMatrixIndexForPoint(std::vector<size_t>{3, 2},
                                                     std::vector<size_t>{7, 3},
                                                     std::vector<size_t>{3, 2}, pre) == 21);
  pre = testObject.test_preComputeMatrixIndexForPoint(std::vector<size_t>{3, 2, 2});
  BOOST_CHECK(testObject.test_getMatrixIndexForPoint(std::vector<size_t>{2, 1, 1},
                                                     std::vector<size_t>{1, 1, 1},
                                                     std::vector<size_t>{3, 2, 2}, pre) == 2);
  BOOST_CHECK(testObject.test_getMatrixIndexForPoint(std::vector<size_t>{2, 2, 1},
                                                     std::vector<size_t>{1, 3, 1},
                                                     std::vector<size_t>{3, 2, 2}, pre) == 11);
  BOOST_CHECK(testObject.test_getMatrixIndexForPoint(std::vector<size_t>{1, 1, 2},
                                                     std::vector<size_t>{1, 1, 1},
                                                     std::vector<size_t>{3, 2, 2}, pre) == 22);
  BOOST_CHECK(testObject.test_getMatrixIndexForPoint(std::vector<size_t>{1, 2, 2},
                                                     std::vector<size_t>{1, 1, 1},
                                                     std::vector<size_t>{3, 2, 2}, pre) == 36);
  BOOST_CHECK(testObject.test_getMatrixIndexForPoint(std::vector<size_t>{3, 2, 2},
                                                     std::vector<size_t>{7, 3, 3},
                                                     std::vector<size_t>{3, 2, 2}, pre) == 63);
}

BOOST_AUTO_TEST_CASE(LhsMatrixPermutationTest) {
  DBMatOfflineTest baseTestObject;
  DBMatOfflineTest desiredTestObject;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.lambda_ = 0;

  sgpp::base::GeneralGridConfiguration baseConfig;
  baseConfig.generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
  baseConfig.type_ = sgpp::base::GridType::Linear;
  baseConfig.dim_ = 3;
  baseConfig.levelVector_ = std::vector<size_t>{3, 2, 2};

  sgpp::base::GeneralGridConfiguration desiredConfig;
  desiredConfig.generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
  desiredConfig.type_ = sgpp::base::GridType::Linear;
  desiredConfig.dim_ = 5;
  desiredConfig.levelVector_ = std::vector<size_t>{3, 2, 1, 1, 2};

  sgpp::datadriven::GridFactory gridFactory;
  std::set<std::set<size_t>> interactions = std::set<std::set<size_t>>();
  std::unique_ptr<sgpp::base::Grid> baseGrid{gridFactory.createGrid(baseConfig, interactions)};
  std::unique_ptr<sgpp::base::Grid> desiredGrid{
      gridFactory.createGrid(desiredConfig, interactions)};

  baseTestObject.buildMatrix(baseGrid.get(), regularizationConfig);
  desiredTestObject.buildMatrix(desiredGrid.get(), regularizationConfig);

  baseTestObject.permuteLhsMatrix(baseConfig, desiredConfig);

  for (size_t i = 0; i < baseTestObject.getLhsMatrix_ONLY_FOR_TESTING().getNrows(); i++) {
    for (size_t j = 0; j < baseTestObject.getLhsMatrix_ONLY_FOR_TESTING().getNcols(); j++) {
      BOOST_CHECK(std::abs(baseTestObject.getLhsMatrix_ONLY_FOR_TESTING().get(i, j) -
                           desiredTestObject.getLhsMatrix_ONLY_FOR_TESTING().get(i, j)) < 0.001);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
