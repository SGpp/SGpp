/*
 * multiEvalPerformance.cpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */
#ifdef __AVX__

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <random>
#include <fstream>
#include <iostream>

#include <zlib.h>


#include "test_datadrivenCommon.hpp"

#include <sgpp/globaldef.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/tools/ConfigurationParameters.hpp>

BOOST_AUTO_TEST_SUITE(TestStreamingMult)

BOOST_AUTO_TEST_CASE(Simple) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E-24), std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 1E-21) };

    uint32_t level = 5;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::DEFAULT);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(SGPP::base::GridType::Linear, std::get<0>(fileNameError), level, configuration);
//        BOOST_TEST_MESSAGE(std::get<0>(fileNameError));
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
