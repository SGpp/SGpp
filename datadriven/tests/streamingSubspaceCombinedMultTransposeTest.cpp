// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#if USE_OCL == 1
#ifdef __AVX__

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <zlib.h>

#include <random>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "test_datadrivenCommon.hpp"
#include "sgpp/globaldef.hpp"
#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "sgpp/base/tools/ConfigurationParameters.hpp"

BOOST_AUTO_TEST_SUITE(TestStreamingSubspaceCombinedMultTranspose)

BOOST_AUTO_TEST_CASE(Simple) {
  std::vector<std::tuple<std::string, double> > fileNamesError = {
      std::tuple<std::string, double>("datadriven/tests/data/friedman_4d.arff.gz", 1E-13),
      std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 1E-17)};

  uint32_t level = 5;

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::SUBSPACELINEAR,
      sgpp::datadriven::OperationMultipleEvalSubType::COMBINED);

  for (std::tuple<std::string, double> fileNameError : fileNamesError) {
    double mse = compareToReferenceTranspose(sgpp::base::GridType::Linear,
                                             std::get<0>(fileNameError), level, configuration);
    //        BOOST_TEST_MESSAGE(std::get<0>(fileNameError));
    BOOST_CHECK(mse < std::get<1>(fileNameError));
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
#endif
