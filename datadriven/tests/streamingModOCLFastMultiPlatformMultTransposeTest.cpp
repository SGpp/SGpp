/*
 * multiEvalPerformance.cpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE StreamingModOCLFastMultiPlatformMultTranspose
#include <boost/test/unit_test.hpp>

#include <random>
#include <fstream>
#include <iostream>

#include <zlib.h>

#include <sgpp/globaldef.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/tools/ConfigurationParameters.hpp>
#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>

std::string uncompressFile(std::string fileName) {

    gzFile inFileZ = gzopen(fileName.c_str(), "rb");
    if (inFileZ == NULL) {
        std::cout << "Error: Failed to gzopen file " << fileName << std::endl;
        exit(0);
    }
    unsigned char unzipBuffer[8192];
    unsigned int unzippedBytes;
    std::vector<unsigned char> unzippedData;
    while (true) {
        unzippedBytes = gzread(inFileZ, unzipBuffer, 8192);
        if (unzippedBytes > 0) {
            for (size_t i = 0; i < unzippedBytes; i++) {
                unzippedData.push_back(unzipBuffer[i]);
            }
        } else {
            break;
        }
    }
    gzclose(inFileZ);

    std::stringstream convert;
    for (size_t i = 0; i < unzippedData.size(); i++) {
        convert << unzippedData[i];
    }
    return convert.str();
}

void doAllRefinements(SGPP::base::AdpativityConfiguration& adaptConfig,
SGPP::base::Grid& grid, SGPP::base::GridGenerator& gridGen) {

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(1, 100);

    SGPP::base::DataVector alphaRefine(grid.getSize());

    for (size_t i = 0; i < alphaRefine.getSize(); i++) {
        alphaRefine[i] = dist(mt);
    }

    for (size_t i = 0; i < adaptConfig.numRefinements_; i++) {
        SGPP::base::SurplusRefinementFunctor* myRefineFunc = new SGPP::base::SurplusRefinementFunctor(&alphaRefine,
                adaptConfig.noPoints_, adaptConfig.threshold_);
        gridGen.refine(myRefineFunc);
        size_t oldSize = alphaRefine.getSize();
        alphaRefine.resize(grid.getSize());

        for (size_t j = oldSize; j < alphaRefine.getSize(); j++) {
            alphaRefine[j] = dist(mt);
        }

        delete myRefineFunc;
    }
}

double compareToReference(std::string fileName, size_t level, SGPP::base::AdpativityConfiguration adaptConfig,
SGPP::datadriven::OperationMultipleEvalConfiguration configuration) {

    std::string content = uncompressFile(fileName);
//
    SGPP::datadriven::ARFFTools arffTools;
    SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

//    SGPP::datadriven::Dataset dataset = arffTools.readARFF("friedman_4d.arff");

    SGPP::base::DataMatrix* trainingData = dataset.getTrainingData();

    size_t dim = dataset.getDimension();
    SGPP::base::Grid* grid = SGPP::base::Grid::createModLinearGrid(dim);
    SGPP::base::GridStorage* gridStorage = grid->getStorage();
    BOOST_TEST_MESSAGE("dimensionality:        " << gridStorage->dim());

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

    BOOST_TEST_MESSAGE("creating operation with unrefined grid");
    SGPP::base::OperationMultipleEval* eval =
    SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData, configuration);

    doAllRefinements(adaptConfig, *grid, *gridGen);

    BOOST_TEST_MESSAGE("number of grid points after refinement: " << gridStorage->size());
    BOOST_TEST_MESSAGE("grid set up");

    SGPP::base::DataVector alphaResult(gridStorage->size());

    BOOST_TEST_MESSAGE("preparing operation for refined grid");
    eval->prepare();

    BOOST_TEST_MESSAGE("calculating result");

    eval->multTranspose(dataSizeVector, alphaResult);

    BOOST_TEST_MESSAGE("calculating comparison values...");

    SGPP::base::OperationMultipleEval* evalCompare =
    SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData);

    SGPP::base::DataVector alphaResultCompare(gridStorage->size());
    alphaResultCompare.setAll(0.0);

    evalCompare->multTranspose(dataSizeVector, alphaResultCompare);

    double mse = 0.0;

    for (size_t i = 0; i < alphaResultCompare.getSize(); i++) {
//        std::cout << "mine: " << alphaResult[i] << " ref: " << alphaResultCompare[i] << std::endl;
        mse += (alphaResult[i] - alphaResultCompare[i]) * (alphaResult[i] - alphaResultCompare[i]);
    }

    mse = mse / static_cast<double>(alphaResultCompare.getSize());
    BOOST_TEST_MESSAGE("mse: " << mse);
    return mse;
}

BOOST_AUTO_TEST_CASE(Simple) {

    //  std::string fileName = "friedman2_90000.arff";
    //    std::string fileName = "debugging.arff";
    std::vector<std::string> fileNames = { "friedman_4d.arff.gz", "friedman_10d.arff.gz" };

    uint32_t level = 4;
    //  uint32_t level = 3;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
    configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["OCL_MANAGER_VERBOSE"] = "false";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    double mse;

    parameters["KERNEL_USE_LOCAL_MEMORY"] = "true";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "1";
    parameters["KERNEL_TRANS_GRID_BLOCK_SIZE"] = "1";
    parameters["KERNEL_TRANS_DATA_BLOCK_SIZE"] = "1";
    parameters["KERNEL_MAX_DIM_UNROLL"] = "1";
    parameters["KERNEL_STORE_DATA"] = "array";
    parameters["PLATFORM"] = "NVIDIA CUDA";
    parameters["SELECT_SPECIFIC_DEVICE"] = "0";
    for (std::string fileName : fileNames) {
        mse = compareToReference(fileName, level, adaptConfig, configuration);
        BOOST_CHECK(mse < 10E-14);
    }
}

BOOST_AUTO_TEST_CASE(Blocking) {

    //  std::string fileName = "friedman2_90000.arff";
    //    std::string fileName = "debugging.arff";
    std::vector<std::string> fileNames = { "friedman_4d.arff.gz", "friedman_10d.arff.gz" };

    uint32_t level = 4;
    //  uint32_t level = 3;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
    configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["OCL_MANAGER_VERBOSE"] = "false";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    double mse;

    parameters["KERNEL_USE_LOCAL_MEMORY"] = "true";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "2";
    parameters["KERNEL_TRANS_GRID_BLOCK_SIZE"] = "2";
    parameters["KERNEL_TRANS_DATA_BLOCK_SIZE"] = "2";
    parameters["KERNEL_MAX_DIM_UNROLL"] = "1";
    parameters["KERNEL_STORE_DATA"] = "array";
    parameters["PLATFORM"] = "NVIDIA CUDA";
    parameters["SELECT_SPECIFIC_DEVICE"] = "0";
    for (std::string fileName : fileNames) {
        mse = compareToReference(fileName, level, adaptConfig, configuration);
        BOOST_CHECK(mse < 10E-14);
    }
}

BOOST_AUTO_TEST_CASE(MultiDevice) {

    //  std::string fileName = "friedman2_90000.arff";
    //    std::string fileName = "debugging.arff";
    std::vector<std::string> fileNames = { "friedman_4d.arff.gz", "friedman_10d.arff.gz" };

    uint32_t level = 4;
    //  uint32_t level = 3;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
    configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["OCL_MANAGER_VERBOSE"] = "false";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    double mse;

    //    OCL_MANAGER_VERBOSE=true
    //    LINEAR_LOAD_BALANCING_VERBOSE=true
    //    LOCAL_SIZE=128
    //    KERNEL_USE_LOCAL_MEMORY=false
    //    KERNEL_DATA_BLOCKING_SIZE=2
    //    REUSE_SOURCE=false
    //    ENABLE_OPTIMIZATIONS=true
    //    SHOW_BUILD_LOG=true
    //    KERNEL_TRANS_DATA_BLOCK_SIZE=2
    //    KERNEL_TRANS_GRID_BLOCK_SIZE=2
    //    KERNEL_MAX_DIM_UNROLL=1
    //    WRITE_SOURCE=true
    //    KERNEL_STORE_DATA=array
    //    PLATFORM=Intel(R) OpenCL
    //    SELECT_SPECIFIC_DEVICE=DISABLED
    //    INTERNAL_PRECISION=double
    //    KERNEL_VERBOSE=true

    //    parameters["OCL_MANAGER_VERBOSE"] = "true";
    //    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "true";
    //    parameters["LOCAL_SIZE"] = "128";
    //    parameters["KERNEL_USE_LOCAL_MEMORY"] = "false";
    //    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "2";
    //    parameters["REUSE_SOURCE"] = "false";
    //    parameters["ENABLE_OPTIMIZATIONS"] = "true";
    //    parameters["SHOW_BUILD_LOG"] = "true";
    //    parameters["KERNEL_TRANS_DATA_BLOCK_SIZE"] = "2";
    //    parameters["KERNEL_TRANS_GRID_BLOCK_SIZE"] = "2";
    //    parameters["KERNEL_MAX_DIM_UNROLL"] = "1";
    //    parameters["WRITE_SOURCE"] = "true";
    //    parameters["KERNEL_STORE_DATA"] = "array";
    //    parameters["PLATFORM"] = "Intel(R) OpenCL";
    //    parameters["SELECT_SPECIFIC_DEVICE"] = "DISABLED";
    //    parameters["INTERNAL_PRECISION"] = "double";
    //    parameters["KERNEL_VERBOSE"] = "true";

    parameters["KERNEL_USE_LOCAL_MEMORY"] = "false";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "2";
    parameters["KERNEL_TRANS_GRID_BLOCK_SIZE"] = "2";
    parameters["KERNEL_TRANS_DATA_BLOCK_SIZE"] = "2";
    parameters["KERNEL_MAX_DIM_UNROLL"] = "1";
    parameters["KERNEL_STORE_DATA"] = "array";
    parameters["PLATFORM"] = "NVIDIA CUDA";
    parameters["SELECT_SPECIFIC_DEVICE"] = "DISABLED";
    for (std::string fileName : fileNames) {
        mse = compareToReference(fileName, level, adaptConfig, configuration);
        BOOST_CHECK(mse < 10E-14);
    }
}

BOOST_AUTO_TEST_CASE(MultiPlatform) {

    //  std::string fileName = "friedman2_90000.arff";
    //    std::string fileName = "debugging.arff";
    std::vector<std::string> fileNames = { "friedman_4d.arff.gz", "friedman_10d.arff.gz" };

    uint32_t level = 4;
    //  uint32_t level = 3;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
    configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["OCL_MANAGER_VERBOSE"] = "false";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    double mse;

    parameters["KERNEL_USE_LOCAL_MEMORY"] = "true";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "2";
    parameters["KERNEL_TRANS_GRID_BLOCK_SIZE"] = "2";
    parameters["KERNEL_TRANS_DATA_BLOCK_SIZE"] = "2";
    parameters["KERNEL_MAX_DIM_UNROLL"] = "1";
    parameters["KERNEL_STORE_DATA"] = "array";
    parameters["PLATFORM"] = "all";
    parameters["SELECT_SPECIFIC_DEVICE"] = "DISABLED";
    for (std::string fileName : fileNames) {
        mse = compareToReference(fileName, level, adaptConfig, configuration);
        BOOST_CHECK(mse < 10E-14);
    }
}

BOOST_AUTO_TEST_CASE(SimpleSinglePrecision) {

//    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
//            "friedman_4d.arff.gz", 10E-7), std::tuple<std::string, double>("friedman_10d.arff.gz", 10E-2) };

//    std::vector<std::tuple<std::string, double> > fileNamesError = { };

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "friedman_4d.arff.gz", 2.0E+04), std::tuple<std::string, double>("friedman_10d.arff.gz", 1.0E+5) };

    uint32_t level = 4;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
    configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["OCL_MANAGER_VERBOSE"] = "false";
    parameters["KERNEL_VERBOSE"] = "false";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    parameters["INTERNAL_PRECISION"] = "float";
    double mse;

    parameters["KERNEL_USE_LOCAL_MEMORY"] = "false";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "1";
    parameters["KERNEL_TRANS_GRID_BLOCK_SIZE"] = "1";
    parameters["KERNEL_TRANS_DATA_BLOCK_SIZE"] = "1";
    parameters["KERNEL_MAX_DIM_UNROLL"] = "1";
    parameters["KERNEL_STORE_DATA"] = "array";
    parameters["PLATFORM"] = "NVIDIA CUDA";
    parameters["WRITE_SOURCE"] = "true";
    parameters["SELECT_SPECIFIC_DEVICE"] = "0";
    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(BlockingSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "friedman_4d.arff.gz", 2E+4), std::tuple<std::string, double>("friedman_10d.arff.gz", 1E+4) };

    uint32_t level = 4;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
    configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["OCL_MANAGER_VERBOSE"] = "false";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    parameters["INTERNAL_PRECISION"] = "float";
    double mse;

    parameters["KERNEL_USE_LOCAL_MEMORY"] = "true";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "2";
    parameters["KERNEL_TRANS_GRID_BLOCK_SIZE"] = "2";
    parameters["KERNEL_TRANS_DATA_BLOCK_SIZE"] = "2";
    parameters["KERNEL_MAX_DIM_UNROLL"] = "1";
    parameters["KERNEL_STORE_DATA"] = "array";
    parameters["PLATFORM"] = "NVIDIA CUDA";
    parameters["SELECT_SPECIFIC_DEVICE"] = "0";
    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiDeviceSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "friedman_4d.arff.gz", 2E+4), std::tuple<std::string, double>("friedman_10d.arff.gz", 1E+4) };

    uint32_t level = 4;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
    configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["OCL_MANAGER_VERBOSE"] = "false";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    parameters["INTERNAL_PRECISION"] = "float";
    double mse;

    parameters["KERNEL_USE_LOCAL_MEMORY"] = "true";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "2";
    parameters["KERNEL_TRANS_GRID_BLOCK_SIZE"] = "2";
    parameters["KERNEL_TRANS_DATA_BLOCK_SIZE"] = "2";
    parameters["KERNEL_MAX_DIM_UNROLL"] = "1";
    parameters["KERNEL_STORE_DATA"] = "array";
    parameters["PLATFORM"] = "NVIDIA CUDA";
    parameters["SELECT_SPECIFIC_DEVICE"] = "DISABLED";
    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(MultiPlatformSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "friedman_4d.arff.gz", 2E+4), std::tuple<std::string, double>("friedman_10d.arff.gz", 1E+4) };

    uint32_t level = 4;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
    configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["OCL_MANAGER_VERBOSE"] = "false";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    parameters["INTERNAL_PRECISION"] = "float";
    double mse;

    parameters["KERNEL_USE_LOCAL_MEMORY"] = "true";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "2";
    parameters["KERNEL_TRANS_GRID_BLOCK_SIZE"] = "2";
    parameters["KERNEL_TRANS_DATA_BLOCK_SIZE"] = "2";
    parameters["KERNEL_MAX_DIM_UNROLL"] = "1";
    parameters["KERNEL_STORE_DATA"] = "array";
    parameters["PLATFORM"] = "all";
    parameters["SELECT_SPECIFIC_DEVICE"] = "DISABLED";
    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}
