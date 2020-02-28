// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflinePermutable.hpp>

#include <vector>

BOOST_AUTO_TEST_SUITE(PermutationUtilTest)

BOOST_AUTO_TEST_CASE(IsPermutationTest) {
  std::vector<size_t> baseVec{2, 3, 2};
  std::vector<size_t> permVec{3, 2, 2};
  BOOST_CHECK(sgpp::datadriven::PermutationUtil::isPermutation(baseVec, permVec));

  baseVec = std::vector<size_t>{2, 3, 4, 5, 6};
  permVec = std::vector<size_t>{6, 5, 4, 3, 2};
  BOOST_CHECK(sgpp::datadriven::PermutationUtil::isPermutation(baseVec, permVec));

  baseVec = std::vector<size_t>{2};
  permVec = std::vector<size_t>{2};
  BOOST_CHECK(sgpp::datadriven::PermutationUtil::isPermutation(baseVec, permVec));

  baseVec = std::vector<size_t>{2, 3, 4};
  permVec = std::vector<size_t>{2, 3, 4};
  BOOST_CHECK(sgpp::datadriven::PermutationUtil::isPermutation(baseVec, permVec));

  baseVec = std::vector<size_t>{2, 3, 3};
  permVec = std::vector<size_t>{3, 3, 3};
  BOOST_CHECK(!sgpp::datadriven::PermutationUtil::isPermutation(baseVec, permVec));

  baseVec = std::vector<size_t>{3, 4, 4};
  permVec = std::vector<size_t>{3, 4, 4, 1};
  BOOST_CHECK(!sgpp::datadriven::PermutationUtil::isPermutation(baseVec, permVec));
}

BOOST_AUTO_TEST_CASE(DeleteOnesFromLevelVecTest) {
  std::vector<size_t> vec{1, 2, 1, 3, 1, 4, 1, 5, 1};
  std::vector<size_t> check{2, 3, 4, 5};
  BOOST_CHECK(sgpp::datadriven::PermutationUtil::deleteOnesFromLevelVec(vec) == check);

  vec = std::vector<size_t>{2, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  check = std::vector<size_t>{2, 3, 3};
  BOOST_CHECK(sgpp::datadriven::PermutationUtil::deleteOnesFromLevelVec(vec) == check);

  vec = std::vector<size_t>{1, 1, 1, 1, 2, 3, 3};
  check = std::vector<size_t>{2, 3, 3};
  BOOST_CHECK(sgpp::datadriven::PermutationUtil::deleteOnesFromLevelVec(vec) == check);
}

BOOST_AUTO_TEST_CASE(NormalizeConfigTest) {
  sgpp::base::GeneralGridConfiguration config;
  config.levelVector_ = {1, 2, 2, 3, 4, 1};
  config.dim_ = 6;
  sgpp::base::GeneralGridConfiguration check =
      sgpp::datadriven::PermutationUtil::getNormalizedConfig(config);

  BOOST_CHECK(check.dim_ == 4 && (check.levelVector_ == std::vector<size_t>{2, 2, 3, 4}));

  config.levelVector_ = {2, 1, 1, 1, 4, 3, 1};
  config.dim_ = 7;
  check =
      sgpp::datadriven::PermutationUtil::getNormalizedConfig(config);

  BOOST_CHECK(check.dim_ == 3 && (check.levelVector_ == std::vector<size_t>{2, 4, 3}));
}

BOOST_AUTO_TEST_SUITE_END()
