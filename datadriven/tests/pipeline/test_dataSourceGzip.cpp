// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef ZLIB
#ifdef USE_GSL

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/ClassificationMinerFactory.hpp>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

using sgpp::datadriven::ClassificationMinerFactory;
using sgpp::datadriven::SparseGridMiner;

BOOST_AUTO_TEST_SUITE(test_gzip_miner)

// Based on the ClassificationMinerFromConfigFile example
BOOST_AUTO_TEST_CASE(classificationMinerGzipTest) {
  const std::string path = "datadriven/tests/pipeline/config_dataSourceGzip.json";
  ClassificationMinerFactory factory;
  auto miner = std::unique_ptr<SparseGridMiner>(factory.buildMiner(path));
  miner->learn(true);
  // This suffices to check that the gzipped dataset could be used. Would error-exit otherwise
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
#endif
