// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef ZLIB
#if USE_OCL == 1
#ifdef __AVX__

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

namespace TestStreamingSubspaceCombinedMultFixture {
struct FilesNamesAndErrorFixture {
  FilesNamesAndErrorFixture() {}
  ~FilesNamesAndErrorFixture() {}

  std::vector<std::tuple<std::string, double>> fileNamesErrorDouble = {
      std::tuple<std::string, double>("datadriven/tests/data/friedman2_4d_10000.arff.gz", 1E-25),
      std::tuple<std::string, double>("datadriven/tests/data/friedman1_10d_2000.arff.gz", 1E-25)};

  uint32_t level = 5;
};
}  // namespace TestStreamingSubspaceCombinedMultFixture

BOOST_FIXTURE_TEST_SUITE(TestStreamingSubspaceCombinedMult,
                         TestStreamingSubspaceCombinedMultFixture::FilesNamesAndErrorFixture)

BOOST_AUTO_TEST_CASE(Simple) {
  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::SUBSPACELINEAR,
      sgpp::datadriven::OperationMultipleEvalSubType::COMBINED);

  compareDatasets(fileNamesErrorDouble, sgpp::base::GridType::Linear, level, configuration);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
#endif
#endif
