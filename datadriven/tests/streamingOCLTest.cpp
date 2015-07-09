/*
 * multiEvalPerformance.cpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE StreamingOCL
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
SGPP::base::Grid& grid, SGPP::base::GridGenerator& gridGen,
SGPP::base::DataVector& alpha) {

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(1, 100);

    for (size_t i = 0; i < adaptConfig.numRefinements_; i++) {
        SGPP::base::SurplusRefinementFunctor* myRefineFunc = new SGPP::base::SurplusRefinementFunctor(&alpha,
                adaptConfig.noPoints_, adaptConfig.threshold_);
        gridGen.refine(myRefineFunc);
        size_t oldSize = alpha.getSize();
        alpha.resize(grid.getSize());

        for (size_t j = oldSize; j < alpha.getSize(); j++) {
            alpha[j] = dist(mt);
        }

        delete myRefineFunc;
    }
}

double compareToReference(std::string fileName, size_t level, SGPP::base::AdpativityConfiguration adaptConfig,
SGPP::datadriven::OperationMultipleEvalConfiguration configuration) {

    std::string content = uncompressFile(fileName);

    SGPP::datadriven::ARFFTools arffTools;
    SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

    SGPP::base::DataMatrix* trainingData = dataset.getTrainingData();

    size_t dim = dataset.getDimension();
    SGPP::base::Grid* grid = SGPP::base::Grid::createLinearGrid(dim);
    SGPP::base::GridStorage* gridStorage = grid->getStorage();
    BOOST_TEST_MESSAGE("dimensionality:        " << gridStorage->dim());

    SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
    gridGen->regular(level);
    BOOST_TEST_MESSAGE("number of grid points: " << gridStorage->size());
    BOOST_TEST_MESSAGE("number of data points: " << dataset.getNumberInstances());

    SGPP::base::DataVector alpha(gridStorage->size());

    for (size_t i = 0; i < alpha.getSize(); i++) {
        //alpha[i] = dist(mt);
        alpha[i] = static_cast<double>(i);
    }

    BOOST_TEST_MESSAGE("creating operation with unrefined grid");
    SGPP::base::OperationMultipleEval* eval =
    SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData, configuration);

    doAllRefinements(adaptConfig, *grid, *gridGen, alpha);

    BOOST_TEST_MESSAGE("number of grid points after refinement: " << gridStorage->size());
    BOOST_TEST_MESSAGE("grid set up");

    SGPP::base::DataVector dataSizeVectorResult(dataset.getNumberInstances());
    dataSizeVectorResult.setAll(0);

    BOOST_TEST_MESSAGE("preparing operation for refined grid");
    eval->prepare();

    BOOST_TEST_MESSAGE("calculating result");

    eval->mult(alpha, dataSizeVectorResult);

    BOOST_TEST_MESSAGE("calculating comparison values...");

    SGPP::base::OperationMultipleEval* evalCompare =
    SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData);

    SGPP::base::DataVector dataSizeVectorResultCompare(dataset.getNumberInstances());
    dataSizeVectorResultCompare.setAll(0.0);

    evalCompare->mult(alpha, dataSizeVectorResultCompare);

    double mse = 0.0;

    for (size_t i = 0; i < dataSizeVectorResultCompare.getSize(); i++) {
//        BOOST_TEST_MESSAGE("mine: " << dataSizeVectorResult[i] << " ref: " << dataSizeVectorResultCompare[i]);
        mse += (dataSizeVectorResult[i] - dataSizeVectorResultCompare[i])
                * (dataSizeVectorResult[i] - dataSizeVectorResultCompare[i]);
    }

    mse = mse / static_cast<double>(dataSizeVectorResultCompare.getSize());
    BOOST_TEST_MESSAGE("mse: " << mse);
    return mse;
}

BOOST_AUTO_TEST_CASE(Simple) {

    //  std::string fileName = "friedman2_90000.arff";
    //    std::string fileName = "debugging.arff";
    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "friedman_4d.arff.gz", 1E-26), std::tuple<std::string, double>("friedman_10d.arff.gz", 1E-24) };

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
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCL;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["OCL_MANAGER_VERBOSE"] = "false";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    double mse;

    parameters["KERNEL_USE_LOCAL_MEMORY"] = "false";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "1";
    parameters["KERNEL_TRANS_GRID_BLOCKING_SIZE"] = "1";
    parameters["KERNEL_STORE_DATA"] = "array";
    parameters["PLATFORM"] = "NVIDIA CUDA";
    parameters["SELECT_SPECIFIC_DEVICE"] = "0";
    parameters["MAX_DEVICES"] = "1";
    parameters["WRITE_SOURCE"] = "true";
    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

/*
BOOST_AUTO_TEST_CASE(Blocking) {

    //  std::string fileName = "friedman2_90000.arff";
//    std::vector<std::string> fileNames = { "debugging.arff.gz" };
    std::vector<std::string> fileNames = { "friedman_4d_small.arff.gz" };
//    std::vector<std::string> fileNames = { "friedman_4d.arff.gz", "friedman_10d.arff.gz" };

    uint32_t level = 1;
    //  uint32_t level = 3;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
    configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCL;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["LOCAL_SIZE"] = "1";
    parameters["OCL_MANAGER_VERBOSE"] = "true";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    parameters["KERNEL_USE_LOCAL_MEMORY"] = "false";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "2";
    parameters["KERNEL_TRANS_GRID_BLOCKING_SIZE"] = "1";
    parameters["KERNEL_STORE_DATA"] = "register";
    parameters["PLATFORM"] = "NVIDIA CUDA";
    parameters["MAX_DEVICES"] = "1";
    parameters["WRITE_SOURCE"] = "true";
    parameters["KERNEL_MAX_DIM_UNROLL"] = "10";
    parameters["SELECT_SPECIFIC_DEVICE"] = "0";
    parameters["REUSE_SOURCE"] = "false";
    for (std::string fileName : fileNames) {
        double mse = compareToReference(fileName, level, adaptConfig, configuration);
        BOOST_CHECK(mse < 10E-14);
    }
}
*/

BOOST_AUTO_TEST_CASE(MultiDevice) {

    //  std::string fileName = "friedman2_90000.arff";
    //    std::string fileName = "debugging.arff";
    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "friedman_4d.arff.gz", 1E-26), std::tuple<std::string, double>("friedman_10d.arff.gz", 1E-24) };

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
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCL;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["OCL_MANAGER_VERBOSE"] = "false";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";

    double mse;

    parameters["KERNEL_USE_LOCAL_MEMORY"] = "false";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "1";
    parameters["KERNEL_TRANS_GRID_BLOCKING_SIZE"] = "1";
    parameters["KERNEL_STORE_DATA"] = "array";
    parameters["PLATFORM"] = "NVIDIA CUDA";
    parameters["MAX_DEVICES"] = "1";
    parameters["SELECT_SPECIFIC_DEVICE"] = "DISABLED";
    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

BOOST_AUTO_TEST_CASE(SimpleSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "friedman_4d.arff.gz", 10E-7), std::tuple<std::string, double>("friedman_10d.arff.gz", 10E-2) };

    uint32_t level = 4;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
    configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCL;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["OCL_MANAGER_VERBOSE"] = "false";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    parameters["INTERNAL_PRECISION"] = "float";
    double mse;

    parameters["KERNEL_USE_LOCAL_MEMORY"] = "false";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "1";
    parameters["KERNEL_TRANS_GRID_BLOCKING_SIZE"] = "1";
    parameters["KERNEL_STORE_DATA"] = "array";
    parameters["PLATFORM"] = "NVIDIA CUDA";
    parameters["MAX_DEVICES"] = "1";
    parameters["SELECT_SPECIFIC_DEVICE"] = "0";
    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}

/*
BOOST_AUTO_TEST_CASE(BlockingSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "friedman_4d.arff.gz", 10E-7), std::tuple<std::string, double>("friedman_10d.arff.gz", 10E-2) };

    uint32_t level = 4;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
    configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCL;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["OCL_MANAGER_VERBOSE"] = "false";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    parameters["INTERNAL_PRECISION"] = "float";
    double mse;

    parameters["KERNEL_USE_LOCAL_MEMORY"] = "false";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "2";
    parameters["KERNEL_TRANS_GRID_BLOCKING_SIZE"] = "1";
    parameters["KERNEL_STORE_DATA"] = "array";
    parameters["PLATFORM"] = "NVIDIA CUDA";
    parameters["MAX_DEVICES"] = "1";
    parameters["SELECT_SPECIFIC_DEVICE"] = "0";
    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}
*/

BOOST_AUTO_TEST_CASE(MultiDeviceSinglePrecision) {

    std::vector<std::tuple<std::string, double> > fileNamesError = { std::tuple<std::string, double>(
            "friedman_4d.arff.gz", 10E-7), std::tuple<std::string, double>("friedman_10d.arff.gz", 10E-2) };

    uint32_t level = 4;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
    configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCL;

    SGPP::base::OCLConfigurationParameters parameters;
    configuration.parameters = (SGPP::base::ConfigurationParameters *) &parameters;
    parameters["OCL_MANAGER_VERBOSE"] = "false";
    parameters["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    parameters["INTERNAL_PRECISION"] = "float";
    double mse;

    parameters["KERNEL_USE_LOCAL_MEMORY"] = "false";
    parameters["KERNEL_DATA_BLOCKING_SIZE"] = "1";
    parameters["KERNEL_TRANS_GRID_BLOCKING_SIZE"] = "1";
    parameters["KERNEL_STORE_DATA"] = "array";
    parameters["PLATFORM"] = "NVIDIA CUDA";
    parameters["MAX_DEVICES"] = "1";
    parameters["SELECT_SPECIFIC_DEVICE"] = "DISABLED";
    for (std::tuple<std::string, double> fileNameError : fileNamesError) {
        mse = compareToReference(std::get<0>(fileNameError), level, adaptConfig, configuration);
        BOOST_CHECK(mse < std::get<1>(fileNameError));
    }
}
