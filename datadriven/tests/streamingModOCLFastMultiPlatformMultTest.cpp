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
#include "../src/sgpp/datadriven/opencl/OCLConfigurationParameters.hpp"

BOOST_AUTO_TEST_SUITE(TestStreamingModOCLFastMultiPlatformMult)

double compareToReference(std::string fileName, size_t level, SGPP::base::AdpativityConfiguration adaptConfig,
SGPP::datadriven::OperationMultipleEvalConfiguration configuration) {

    std::string content = uncompressFile(fileName);

    SGPP::datadriven::ARFFTools arffTools;
    SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

    SGPP::base::DataMatrix* trainingData = dataset.getTrainingData();

    size_t dim = dataset.getDimension();
    SGPP::base::Grid* grid = SGPP::base::Grid::createModLinearGrid(dim);
    SGPP::base::GridStorage* gridStorage = grid->getStorage();
//    BOOST_TEST_MESSAGE("dimensionality:        " << gridStorage->dim());

    SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
    gridGen->regular(level);
//    BOOST_TEST_MESSAGE("number of grid points: " << gridStorage->size());
//    BOOST_TEST_MESSAGE("number of data points: " << dataset.getNumberInstances());

    SGPP::base::DataVector alpha(gridStorage->size());

    for (size_t i = 0; i < alpha.getSize(); i++) {
        //alpha[i] = dist(mt);
        alpha[i] = static_cast<double>(i);
    }

//    BOOST_TEST_MESSAGE("creating operation with unrefined grid");
    SGPP::base::OperationMultipleEval* eval =
    SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData, configuration);

    doRandomRefinements(adaptConfig, *grid, *gridGen, alpha);

//    BOOST_TEST_MESSAGE("number of grid points after refinement: " << gridStorage->size());
//    BOOST_TEST_MESSAGE("grid set up");

    SGPP::base::DataVector dataSizeVectorResult(dataset.getNumberInstances());
    dataSizeVectorResult.setAll(0);

//    BOOST_TEST_MESSAGE("preparing operation for refined grid");
    eval->prepare();

//    BOOST_TEST_MESSAGE("calculating result");

    eval->mult(alpha, dataSizeVectorResult);

//    BOOST_TEST_MESSAGE("calculating comparison values...");

    SGPP::base::OperationMultipleEval* evalCompare =
    SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData);

    SGPP::base::DataVector dataSizeVectorResultCompare(dataset.getNumberInstances());
    dataSizeVectorResultCompare.setAll(0.0);

    evalCompare->mult(alpha, dataSizeVectorResultCompare);

    double mse = 0.0;

    for (size_t i = 0; i < dataSizeVectorResultCompare.getSize(); i++) {
        //std::cout << "mine: " << dataSizeVectorResult[i] << " ref: " << dataSizeVectorResultCompare[i] << std::endl;
        mse += (dataSizeVectorResult[i] - dataSizeVectorResultCompare[i])
                * (dataSizeVectorResult[i] - dataSizeVectorResultCompare[i]);
    }

    mse = mse / static_cast<double>(dataSizeVectorResultCompare.getSize());
//    BOOST_TEST_MESSAGE("mse: " << mse);
    return mse;
}

BOOST_AUTO_TEST_CASE(Simple) {

    //  std::string fileName = "friedman2_90000.arff";
    //    std::string fileName = "debugging.arff";
    std::vector<std::string> fileNames = { "datadriven/tests/data/friedman_4d.arff.gz",
            "datadriven/tests/data/friedman_10d.arff.gz" };

    uint32_t level = 1;
    //  uint32_t level = 3;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "true");
    parameters.set("VERBOSE", "true");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "1");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "1");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "1");
    parameters.set("KERNEL_STORE_DATA", "array");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "0");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::string fileName : fileNames) {
        double mse = compareToReference(fileName, level, adaptConfig, configuration);
        BOOST_CHECK(mse < 10E-14);
    }
}

BOOST_AUTO_TEST_CASE(Blocking) {

    //  std::string fileName = "friedman2_90000.arff";
    //    std::string fileName = "debugging.arff";
    std::vector<std::string> fileNames = { "datadriven/tests/data/friedman_4d.arff.gz",
            "datadriven/tests/data/friedman_10d.arff.gz" };

    uint32_t level = 4;
    //  uint32_t level = 3;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "2");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "2");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "2");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "1");
    parameters.set("KERNEL_STORE_DATA", "array");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "0");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::string fileName : fileNames) {
        double mse = compareToReference(fileName, level, adaptConfig, configuration);
        BOOST_CHECK(mse < 10E-14);
    }
}

BOOST_AUTO_TEST_CASE(MultiDevice) {

    //  std::string fileName = "friedman2_90000.arff";
    //    std::string fileName = "debugging.arff";
    std::vector<std::string> fileNames = { "datadriven/tests/data/friedman_4d.arff.gz",
            "datadriven/tests/data/friedman_10d.arff.gz" };

    uint32_t level = 4;
    //  uint32_t level = 3;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "false");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "2");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "2");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "2");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "1");
    parameters.set("KERNEL_STORE_DATA", "array");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::string fileName : fileNames) {
        double mse = compareToReference(fileName, level, adaptConfig, configuration);
        BOOST_CHECK(mse < 10E-14);
    }
}

BOOST_AUTO_TEST_CASE(MultiPlatform) {

    //  std::string fileName = "friedman2_90000.arff";
    //    std::string fileName = "debugging.arff";
    std::vector<std::string> fileNames = { "datadriven/tests/data/friedman_4d.arff.gz",
            "datadriven/tests/data/friedman_10d.arff.gz" };

    uint32_t level = 4;
    //  uint32_t level = 3;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
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
        double mse = compareToReference(fileName, level, adaptConfig, configuration);
        BOOST_CHECK(mse < 10E-14);
    }
}

BOOST_AUTO_TEST_CASE(SimpleSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 10E-7), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 10E-2) };

    uint32_t level = 4;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "1");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "1");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "1");
    parameters.set("KERNEL_STORE_DATA", "array");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "0");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(BlockingSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 10E-7), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 10E-2) };

    uint32_t level = 4;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "2");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "2");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "2");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "1");
    parameters.set("KERNEL_STORE_DATA", "array");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "0");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiDeviceSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 10E-7), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 10E-2) };

    uint32_t level = 4;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("INTERNAL_PRECISION", "float");
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "2");
    parameters.set("KERNEL_TRANS_GRID_BLOCK_SIZE", "2");
    parameters.set("KERNEL_TRANS_DATA_BLOCK_SIZE", "2");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "1");
    parameters.set("KERNEL_STORE_DATA", "array");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM, parameters);

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiPlatformSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_4d.arff.gz", 10E-7), std::tuple<std::string, double>(
            "datadriven/tests/data/friedman_10d.arff.gz", 10E-2) };

    uint32_t level = 4;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::base::OCLConfigurationParameters parameters;
    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("INTERNAL_PRECISION", "float");
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

    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        double mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
