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

double compareToReference(std::string fileName, size_t level, SGPP::base::AdpativityConfiguration adaptConfig,
SGPP::datadriven::OperationMultipleEvalConfiguration configuration) {

    std::string content = uncompressFile(fileName);

    SGPP::datadriven::ARFFTools arffTools;
    SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

    SGPP::base::DataMatrix* trainingData = dataset.getTrainingData();

    size_t dim = dataset.getDimension();
    SGPP::base::Grid* grid = SGPP::base::Grid::createLinearGrid(dim);
    SGPP::base::GridStorage* gridStorage = grid->getStorage();
//    BOOST_TEST_MESSAGE("dimensionality:        " << gridStorage->dim());

    SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
    gridGen->regular(level);
    BOOST_TEST_MESSAGE("number of grid points: " << gridStorage->size());
    BOOST_TEST_MESSAGE("number of data points: " << dataset.getNumberInstances());

    SGPP::base::DataVector dataSizeVector(dataset.getNumberInstances());

    //Don't use random data! Random data will change the expected MSE
    for (size_t i = 0; i < dataSizeVector.getSize(); i++) {
        //    dataSizeVector[i] = dist(mt);
        dataSizeVector[i] = static_cast<double>(i + 1);
    }

//    BOOST_TEST_MESSAGE("creating operation with unrefined grid");
    SGPP::base::OperationMultipleEval* eval =
    SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData, configuration);

    doRandomRefinements(adaptConfig, *grid, *gridGen);

//    BOOST_TEST_MESSAGE("number of grid points after refinement: " << gridStorage->size());
//    BOOST_TEST_MESSAGE("grid set up");

    SGPP::base::DataVector alphaResult(gridStorage->size());
    alphaResult.setAll(0);

//    BOOST_TEST_MESSAGE("preparing operation for refined grid");
    eval->prepare();

//    BOOST_TEST_MESSAGE("calculating result");

    eval->multTranspose(dataSizeVector, alphaResult);

//    BOOST_TEST_MESSAGE("calculating comparison values...");

    SGPP::base::OperationMultipleEval* evalCompare =
    SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData);

    SGPP::base::DataVector alphaResultCompare(gridStorage->size());
    alphaResultCompare.setAll(0.0);

    evalCompare->multTranspose(dataSizeVector, alphaResultCompare);

    double mse = 0.0;

    for (size_t i = 0; i < alphaResultCompare.getSize(); i++) {
        BOOST_TEST_MESSAGE("i: " << i << " mine: " << alphaResult[i] << " ref: " << alphaResultCompare[i]);
        mse += (alphaResult[i] - alphaResultCompare[i]) * (alphaResult[i] - alphaResultCompare[i]);
    }

    mse = mse / static_cast<double>(alphaResultCompare.getSize());
    BOOST_TEST_MESSAGE("fileName: " << fileName << " mse: " << mse);
    return mse;
}

BOOST_AUTO_TEST_CASE(Simple) {

    //  std::string fileName = "friedman2_90000.arff";
    //    std::string fileName = "debugging.arff";
//    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
//            "datadriven/tests/data/friedman_4d.arff.gz", 1E-13), std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 1E-17) };
    std::vector<std::tuple<std::string, double> > fileNamesError = {  std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 1E-17) };

    uint32_t level = 5;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("VERBOSE", "true");
    parameters.set("ENABLE_OPTIMIZATIONS", "true");
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
        double mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(Blocking) {

//    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
//            "datadriven/tests/data/friedman_4d.arff.gz", 1E-13), std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 1E-17) };
    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 1E-17) };

    uint32_t level = 6;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("ENABLE_OPTIMIZATIONS", "true");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "false");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_TRANS_GRID_BLOCKING_SIZE", "4");
    parameters.set("KERNEL_STORE_DATA", "register");
    parameters.set("PLATFORM", "first");
    parameters.set("MAX_DEVICES", "1");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("SELECT_SPECIFIC_DEVICE", "0");
    parameters.set("REUSE_SOURCE", "false");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiDevice) {

    //  std::string fileName = "friedman2_90000.arff";
    //    std::string fileName = "debugging.arff";
    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E-13), std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 1E-17) };

    uint32_t level = 6;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "true");
    parameters.set("ENABLE_OPTIMIZATIONS", "true");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_TRANS_GRID_BLOCKING_SIZE", "4");
    parameters.set("KERNEL_STORE_DATA", "register");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");
    parameters.set("VERBOSE", "true");


    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(SimpleSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E+4), std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 20.0) };

    uint32_t level = 6;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("ENABLE_OPTIMIZATIONS", "true");
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
        double mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}


BOOST_AUTO_TEST_CASE(BlockingSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E+4), std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 20.0) };

    uint32_t level = 6;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("ENABLE_OPTIMIZATIONS", "true");
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
        double mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiDeviceSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 1E+4), std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 20.0) };

    uint32_t level = 6;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("ENABLE_OPTIMIZATIONS", "true");
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
        double mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
