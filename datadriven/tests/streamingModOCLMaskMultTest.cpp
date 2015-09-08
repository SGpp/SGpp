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

#include <zlib.h>

#include "test_datadrivenCommon.hpp"

#include <sgpp/globaldef.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/tools/ConfigurationParameters.hpp>
#include <sgpp/datadriven/opencl/OCLConfigurationParameters.hpp>

BOOST_AUTO_TEST_SUITE(TestStreamingOCLMaskMult)

SGPP::base::OCLConfigurationParameters getConfigurationDefaults() {
    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("VERBOSE", "false");
    parameters.set("ENABLE_OPTIMIZATIONS", "true");
    parameters.set("PLATFORM", "first");
    parameters.set("MAX_DEVICES", "1");
    parameters.set("SELECT_SPECIFIC_DEVICE", "0");
    return parameters;
}

//double compareToReference(std::string fileName, size_t level,
//SGPP::datadriven::OperationMultipleEvalConfiguration configuration) {
//
//    SGPP::base::AdpativityConfiguration adaptConfig;
//    adaptConfig.maxLevelType_ = false;
//    adaptConfig.noPoints_ = 80;
//    adaptConfig.numRefinements_ = 1;
//    adaptConfig.percent_ = 200.0;
//    adaptConfig.threshold_ = 0.0;
//
//    std::string content = uncompressFile(fileName);
//
//    SGPP::datadriven::ARFFTools arffTools;
//    SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
//
//    SGPP::base::DataMatrix &trainingData = dataset.getTrainingData();
//
//    size_t dim = dataset.getDimension();
//    SGPP::base::Grid* grid = SGPP::base::Grid::createModLinearGrid(dim);
//    SGPP::base::GridStorage* gridStorage = grid->getStorage();
//
//    SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
//    gridGen->regular(level);
//    BOOST_TEST_MESSAGE("number of grid points: " << gridStorage->size());
//    BOOST_TEST_MESSAGE("number of data points: " << dataset.getNumberInstances());
//
//    SGPP::base::DataVector alpha(gridStorage->size());
//
//    for (size_t i = 0; i < alpha.getSize(); i++) {
//        alpha[i] = static_cast<double>(i);
//    }
//
//    SGPP::base::OperationMultipleEval* eval =
//    SGPP::op_factory::createOperationMultipleEval(*grid, trainingData, configuration);
//
//    eval->prepare();
//
//    doRandomRefinements(adaptConfig, *grid, *gridGen, alpha);
//
//    SGPP::base::DataVector dataSizeVectorResult(dataset.getNumberInstances());
//    dataSizeVectorResult.setAll(0);
//
//    eval->prepare();
//
//    eval->mult(alpha, dataSizeVectorResult);
//
//    SGPP::base::OperationMultipleEval* evalCompare =
//    SGPP::op_factory::createOperationMultipleEval(*grid, trainingData);
//
//    SGPP::base::DataVector dataSizeVectorResultCompare(dataset.getNumberInstances());
//    dataSizeVectorResultCompare.setAll(0.0);
//
//    evalCompare->mult(alpha, dataSizeVectorResultCompare);
//
//    double mse = compareVectors(dataSizeVectorResult, dataSizeVectorResultCompare);
//
//    BOOST_TEST_MESSAGE("fileName: " << fileName << " mse: " << mse);
//    return mse;
//}

BOOST_AUTO_TEST_CASE(Simple) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E-22), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 1E-16) };

    uint32_t level = 5;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "false");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "1");


    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMASK, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(SGPP::base::ModLinear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(Local) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E-20), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 1E-13) };

    uint32_t level = 5;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMASK, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(SGPP::base::ModLinear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiDevice) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E-20), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 1E-13) };

    uint32_t level = 5;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("MAX_DEVICES", "0");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMASK, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(SGPP::base::ModLinear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(SimpleSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 10E-3), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 10E5) };

    uint32_t level = 5;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "false");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "1");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMASK, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(SGPP::base::ModLinear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(LocalSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 10E-3), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 10E5) };

    uint32_t level = 5;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMASK, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(SGPP::base::ModLinear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiDeviceSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 10E-3), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 10E5) };

    uint32_t level = 5;

    SGPP::base::OCLConfigurationParameters parameters = getConfigurationDefaults();
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("MAX_DEVICES", "0");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMASK, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(SGPP::base::ModLinear, std::get<0>(fileNameError), level, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
