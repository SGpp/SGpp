// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <zlib.h>
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/base/tools/ConfigurationParameters.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/globaldef.hpp"
#include "test_datadrivenCommon.hpp"

namespace TestStreamingModMaskMultFixture {
struct FilesNamesAndErrorFixture {
  FilesNamesAndErrorFixture() {}
  ~FilesNamesAndErrorFixture() {}

  std::vector<std::tuple<std::string, double>> fileNamesErrorDouble = {
      std::tuple<std::string, double>("datadriven/datasets/friedman/friedman2_4d_10000.arff.gz", 1E-22),
      std::tuple<std::string, double>("datadriven/datasets/friedman/friedman1_10d_2000.arff.gz", 1E-16)};

  uint32_t level = 5;
};
}  // namespace TestStreamingModMaskMultFixture

BOOST_FIXTURE_TEST_SUITE(TestStreamingModMaskMult,
                         TestStreamingModMaskMultFixture::FilesNamesAndErrorFixture)

BOOST_AUTO_TEST_CASE(Simple) {
  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT);

  compareDatasets(fileNamesErrorDouble, sgpp::base::GridType::ModLinear, level, configuration);
}

BOOST_AUTO_TEST_SUITE_END()
