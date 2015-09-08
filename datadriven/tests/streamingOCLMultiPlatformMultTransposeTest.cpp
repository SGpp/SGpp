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
#include <sgpp/datadriven/opencl/OCLConfigurationParameters.hpp>

BOOST_AUTO_TEST_SUITE(TestStreamingOCLMultiPlatformMultTranspose)

SGPP::base::OCLConfigurationParameters getConfigurationDefaults() {
    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("VERBOSE", "false");
    parameters.set("ENABLE_OPTIMIZATIONS", "true");
    return parameters;
}

BOOST_AUTO_TEST_CASE(Simple) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 1E-17) };

    uint32_t level = 5;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_TRANS_GRID_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "1");
    parameters.set("KERNEL_STORE_DATA", "array");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "0");
    parameters.set("MAX_DEVICES", "1");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReferenceTranspose(SGPP::base::Linear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(Blocking) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 1E-17) };

    uint32_t level = 6;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "false");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_TRANS_GRID_BLOCKING_SIZE", "4");
    parameters.set("KERNEL_STORE_DATA", "register");
    parameters.set("PLATFORM", "first");
    parameters.set("MAX_DEVICES", "1");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("SELECT_SPECIFIC_DEVICE", "0");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReferenceTranspose(SGPP::base::Linear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiDevice) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E-13), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 1E-17) };

    uint32_t level = 6;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_TRANS_GRID_BLOCKING_SIZE", "4");
    parameters.set("KERNEL_STORE_DATA", "register");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReferenceTranspose(SGPP::base::Linear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiPlatform) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E-13), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 1E-17) };

    uint32_t level = 6;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_TRANS_GRID_BLOCKING_SIZE", "4");
    parameters.set("KERNEL_STORE_DATA", "register");
    parameters.set("PLATFORM", "all");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReferenceTranspose(SGPP::base::Linear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(SimpleSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E+4), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 20.0) };

    uint32_t level = 6;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_TRANS_GRID_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_STORE_DATA", "array");
    parameters.set("PLATFORM", "first");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("MAX_DEVICES", "1");
    parameters.set("SELECT_SPECIFIC_DEVICE", "0");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReferenceTranspose(SGPP::base::Linear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(BlockingSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E+4), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 20.0) };

    uint32_t level = 6;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "4");
    parameters.set("KERNEL_TRANS_GRID_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_STORE_DATA", "array");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("PLATFORM", "first");
    parameters.set("MAX_DEVICES", "1");
    parameters.set("SELECT_SPECIFIC_DEVICE", "0");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReferenceTranspose(SGPP::base::Linear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiDeviceSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E+4), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 20.0) };

    uint32_t level = 6;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_TRANS_GRID_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_STORE_DATA", "array");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReferenceTranspose(SGPP::base::Linear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiPlatformSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E+4), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 20.0) };

    uint32_t level = 6;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_TRANS_GRID_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_STORE_DATA", "array");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("PLATFORM", "all");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReferenceTranspose(SGPP::base::Linear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
