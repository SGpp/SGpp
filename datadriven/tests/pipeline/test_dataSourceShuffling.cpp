// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorRandom.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorSequential.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorCrossValidation.hpp>
#include <sgpp/datadriven/configuration/CrossvalidationConfiguration.hpp>

#include <vector>

using sgpp::datadriven::DataShufflingFunctor;
using sgpp::datadriven::DataShufflingFunctorSequential;
using sgpp::datadriven::DataShufflingFunctorRandom;
using sgpp::datadriven::DataShufflingFunctorCrossValidation;
using sgpp::datadriven::CrossvalidationConfiguration;

bool testBijectivity(DataShufflingFunctor& shuffling, size_t numSamples) {
  std::vector<bool> hit(numSamples, false);
  for (size_t idx = 0; idx < numSamples; idx++) {
    size_t mappedTo = shuffling(idx, numSamples);
    if (hit[mappedTo]) return false;
    hit[mappedTo] = true;
  }
  return true;
}

bool testOrder(DataShufflingFunctor& shuffling, std::vector<size_t> expectedOrder) {
  for (size_t idx = 0; idx < expectedOrder.size(); idx++) {
    if (shuffling(idx, expectedOrder.size()) != expectedOrder[idx]) return false;
  }
  return true;
}

BOOST_AUTO_TEST_SUITE(testDataSourceShuffling)

BOOST_AUTO_TEST_CASE(TestShufflingSequential) {
  DataShufflingFunctorSequential shuffling;
  BOOST_CHECK(testBijectivity(shuffling, 5000));
}

BOOST_AUTO_TEST_CASE(TestShufflingRandom) {
  DataShufflingFunctorRandom shuffling{-1};
  BOOST_CHECK(testBijectivity(shuffling, 5000));
}

BOOST_AUTO_TEST_CASE(TestShufflingCrossValidation) {
  DataShufflingFunctorSequential seqShuffling;
  CrossvalidationConfiguration crossValidationConfig;
  crossValidationConfig.kfold_ = 3;
  DataShufflingFunctorCrossValidation cvShuffling(crossValidationConfig, &seqShuffling);

  // Test for fold = 0
  cvShuffling.setFold(0);
  std::vector<size_t> expectedOrder{0, 1, 2, 3, 4, 5, 6};
  BOOST_CHECK(testBijectivity(cvShuffling, 5000));
  BOOST_CHECK(testOrder(cvShuffling, expectedOrder));

  // Test for fold = 1
  cvShuffling.setFold(1);
  expectedOrder = std::vector<size_t>{2, 3, 0, 1, 4, 5, 6};
  BOOST_CHECK(testBijectivity(cvShuffling, 5000));
  BOOST_CHECK(testOrder(cvShuffling, expectedOrder));

  // Test for fold = 2
  cvShuffling.setFold(2);
  expectedOrder = std::vector<size_t>{4, 5, 6, 0, 1, 2, 3};
  BOOST_CHECK(testBijectivity(cvShuffling, 5000));
  BOOST_CHECK(testOrder(cvShuffling, expectedOrder));
}

BOOST_AUTO_TEST_SUITE_END()





