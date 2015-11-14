/*
 * multiEvalPerformance.cpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#if USE_OCL==1

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <random>
#include <fstream>
#include <iostream>

#include "test_datadrivenCommon.hpp"

#include <sgpp/globaldef.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/tools/ConfigurationParameters.hpp>
#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>

BOOST_AUTO_TEST_SUITE(TestStreamingModOCLFastMultiPlatformMult)

SGPP::base::OCLConfigurationParameters getConfigurationDefaults() {
    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("VERBOSE", "false");
    parameters.set("ENABLE_OPTIMIZATIONS", "true");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "0");
    return parameters;
}

BOOST_AUTO_TEST_CASE(Simple) {

//    std::vector<std::string> fileNames = { "datadriven/tests/data/friedman_4d.arff.gz",
//            "datadriven/tests/data/friedman_10d.arff.gz" };

    std::vector<std::string> fileNames = { "datadriven/tests/data/friedman_4d.arff.gz" };

    uint32_t level = 4;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "1");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "1");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "1");
    parameters.set("KERNEL_STORE_DATA", "array");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::string fileName : fileNames) {
        double mse = compareToReference(SGPP::base::GridType::ModLinear, fileName, level, configuration);
        BOOST_CHECK(mse < 10E-14);
    }
}

BOOST_AUTO_TEST_CASE(Blocking) {

    std::vector<std::string> fileNames = { "datadriven/tests/data/friedman_4d.arff.gz",
            "datadriven/tests/data/friedman_10d.arff.gz" };

    uint32_t level = 4;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "false");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "2");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "2");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "2");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("KERNEL_STORE_DATA", "register");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::string fileName : fileNames) {
        double mse = compareToReference(SGPP::base::GridType::ModLinear, fileName, level, configuration);
        BOOST_CHECK(mse < 10E-14);
    }
}

BOOST_AUTO_TEST_CASE(MultiDevice) {

    std::vector<std::string> fileNames = { "datadriven/tests/data/friedman_4d.arff.gz",
            "datadriven/tests/data/friedman_10d.arff.gz" };

    uint32_t level = 4;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "false");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "2");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "2");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "2");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("KERNEL_STORE_DATA", "register");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::string fileName : fileNames) {
        double mse = compareToReference(SGPP::base::GridType::ModLinear, fileName, level, configuration);
        BOOST_CHECK(mse < 10E-14);
    }
}

BOOST_AUTO_TEST_CASE(MultiPlatform) {

    std::vector<std::string> fileNames = { "datadriven/tests/data/friedman_4d.arff.gz",
            "datadriven/tests/data/friedman_10d.arff.gz" };

    uint32_t level = 4;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "2");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "2");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "2");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "1");
    parameters.set("KERNEL_STORE_DATA", "array");
    parameters.set("PLATFORM", "all");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::string fileName : fileNames) {
        double mse = compareToReference(SGPP::base::GridType::ModLinear, fileName, level, configuration);
        BOOST_CHECK(mse < 10E-14);
    }
}

BOOST_AUTO_TEST_CASE(SimpleSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 10E-7), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 10E-2) };

    uint32_t level = 4;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "1");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "1");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("KERNEL_STORE_DATA", "register");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(SGPP::base::GridType::ModLinear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(BlockingSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 10E-7), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 10E-2) };

    uint32_t level = 4;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "2");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "2");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "2");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("KERNEL_STORE_DATA", "register");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(SGPP::base::GridType::ModLinear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiDeviceSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 10E-7), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 10E-2) };

    uint32_t level = 4;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "2");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "2");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "2");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("KERNEL_STORE_DATA", "register");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(SGPP::base::GridType::ModLinear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiPlatformSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 10E-7), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 10E-2) };

    uint32_t level = 4;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "2");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "2");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "2");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("KERNEL_STORE_DATA", "register");
    parameters.set("PLATFORM", "all");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(SGPP::base::GridType::ModLinear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
