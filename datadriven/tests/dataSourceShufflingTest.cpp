/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 * DBMatOfflineDatabaseTest.cpp
 *
 * dataSourceShufflingTest.cpp
 *
 *  Created on: Jul 24, 2018
 *      Author: dominik
 */

#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorRandom.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorSequential.hpp>

#include <vector>

using sgpp::datadriven::DataShufflingFunctor;
using sgpp::datadriven::DataShufflingFunctorSequential;
using sgpp::datadriven::DataShufflingFunctorRandom;

bool testBijectivity(DataShufflingFunctor& shuffling, size_t numSamples) {
  std::vector<bool> hit(numSamples, false);
  for (size_t idx = 0; idx < numSamples; idx++) {
    size_t mappedTo = shuffling(idx, numSamples);
    if (hit[mappedTo]) return false;
    hit[mappedTo] = true;
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

BOOST_AUTO_TEST_SUITE_END()





