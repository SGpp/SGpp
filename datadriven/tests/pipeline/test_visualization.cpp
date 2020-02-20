// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#define BOOST_TEST_DYN_LINK
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/ClassificationMinerFactory.hpp>

#include <string>

using sgpp::datadriven::ClassificationMinerFactory;
using sgpp::datadriven::SparseGridMiner;

BOOST_AUTO_TEST_SUITE(testVisualization)

BOOST_AUTO_TEST_CASE(visualization) {
  std::string configFile = "datadriven/tests/pipeline/config_visualization.json";
  ClassificationMinerFactory factory;
  SparseGridMiner *miner = factory.buildMiner(configFile);
  miner->learn(false);
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* USE_GSL */
