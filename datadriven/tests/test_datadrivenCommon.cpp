// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef ZLIB

#include "test_datadrivenCommon.hpp"

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalScalapack/OperationMultipleEvalDistributed.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <zlib.h>
#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::datadriven::BlacsProcessGrid;
using sgpp::datadriven::DataVectorDistributed;
using sgpp::datadriven::OperationMultipleEvalDistributed;

#if USE_OCL == 1
#include "sgpp/base/opencl/OCLManagerMultiPlatform.hpp"
using sgpp::base::OCLManagerMultiPlatform;
using sgpp::base::OCLOperationConfiguration;
#endif

std::string uncompressFile(std::string fileName) {
  gzFile inFileZ = gzopen(fileName.c_str(), "rb");

  if (inFileZ == nullptr) {
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

DataMatrix* readReferenceMatrix(sgpp::base::GridStorage& storage, std::string fileName) {
  std::string content = uncompressFile(fileName);

  std::stringstream contentStream;
  contentStream << content;
  std::string line;

  DataMatrix* m = new DataMatrix(0, storage.getSize());

  size_t currentRow = 0;

  while (!contentStream.eof()) {
    std::getline(contentStream, line);

    // for lines that only contain a newline
    if (line.size() == 0) {
      break;
    }

    m->appendRow();

    size_t curPos = 0;
    size_t curFind = 0;
    std::string curValue;
    double floatValue;

    for (size_t i = 0; i < storage.getSize(); i++) {
      curFind = line.find_first_of(" \t", curPos);
      curValue = line.substr(curPos, curFind - curPos);

      floatValue = boost::lexical_cast<double>(curValue);
      m->set(currentRow, i, floatValue);
      curPos = curFind + 1;
    }

    currentRow += 1;
  }

  return m;
}

void doRandomRefinements(sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                         sgpp::base::Grid& grid, sgpp::base::GridGenerator& gridGen,
                         sgpp::base::DataVector& alpha) {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  for (size_t i = 0; i < adaptivityConfig.numRefinements_; i++) {
    sgpp::base::SurplusRefinementFunctor myRefineFunc(alpha, adaptivityConfig.numRefinementPoints_,
                                                      adaptivityConfig.refinementThreshold_);
    gridGen.refine(myRefineFunc);
    size_t oldSize = alpha.getSize();
    alpha.resize(grid.getSize());

    for (size_t j = oldSize; j < alpha.getSize(); j++) {
      alpha[j] = dist(mt);
    }
  }
}

void doRandomRefinements(sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                         sgpp::base::Grid& grid, sgpp::base::GridGenerator& gridGen) {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  sgpp::base::DataVector alphaRefine(grid.getSize());

  for (size_t i = 0; i < alphaRefine.getSize(); i++) {
    alphaRefine[i] = dist(mt);
  }

  for (size_t i = 0; i < adaptivityConfig.numRefinements_; i++) {
    sgpp::base::SurplusRefinementFunctor myRefineFunc(
        alphaRefine, adaptivityConfig.numRefinementPoints_, adaptivityConfig.refinementThreshold_);
    gridGen.refine(myRefineFunc);
    size_t oldSize = alphaRefine.getSize();
    alphaRefine.resize(grid.getSize());

    for (size_t j = oldSize; j < alphaRefine.getSize(); j++) {
      alphaRefine[j] = dist(mt);
    }
  }
}

double compareVectors(sgpp::base::DataVector& results, sgpp::base::DataVector& resultsCompare) {
  double mse = 0.0;

  bool anyDifferentValue = false;
  double largestDifference = 0.0;
  double value = 0.0;
  double valueReference = 0.0;

  for (size_t i = 0; i < resultsCompare.getSize(); i++) {
    double diff = (results[i] - resultsCompare[i]);

    if (diff > largestDifference) {
      anyDifferentValue = true;
      largestDifference = diff;
      value = results[i];
      valueReference = resultsCompare[i];
    }

    mse += diff * diff;
  }

  if (anyDifferentValue) {
    BOOST_TEST_MESSAGE("largestDifference: " << largestDifference << " value: " << value
                                             << " valueReference: " << valueReference);
  } else {
    BOOST_TEST_MESSAGE("every value matched exactly");
  }

  mse = mse / static_cast<double>(resultsCompare.getSize());
  return mse;
}

void compareDatasets(const std::vector<std::tuple<std::string, double>>& fileNamesError,
                     sgpp::base::GridType gridType, size_t level,
                     sgpp::datadriven::OperationMultipleEvalConfiguration configuration) {
  for (std::tuple<std::string, double> fileNameError : fileNamesError) {
    double mse = compareToReference(gridType, std::get<0>(fileNameError), level, configuration);
    BOOST_CHECK(mse < std::get<1>(fileNameError));
    std::cout << "expected error: " << std::get<1>(fileNameError) << ", observed error:" << mse
              << std::endl;
  }
}

double compareToReference(sgpp::base::GridType gridType, const std::string& fileName, size_t level,
                          sgpp::datadriven::OperationMultipleEvalConfiguration configuration) {
  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.numRefinementPoints_ = 80;
  adaptivityConfig.numRefinements_ = 1;
  adaptivityConfig.percent_ = 200.0;
  adaptivityConfig.refinementThreshold_ = 0.0;

  std::string content = uncompressFile(fileName);

  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

  sgpp::base::DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::shared_ptr<sgpp::base::Grid> grid;

  if (gridType == sgpp::base::GridType::Linear) {
    grid = std::shared_ptr<sgpp::base::Grid>(sgpp::base::Grid::createLinearGrid(dim));
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    grid = std::shared_ptr<sgpp::base::Grid>(sgpp::base::Grid::createModLinearGrid(dim));
  }

  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(level);

  sgpp::base::DataVector alpha(gridStorage.getSize());

  for (size_t i = 0; i < alpha.getSize(); i++) {
    alpha[i] = static_cast<double>(i);
  }

  auto eval = std::shared_ptr<sgpp::base::OperationMultipleEval>(
      sgpp::op_factory::createOperationMultipleEval(*grid, trainingData, configuration));

  eval->prepare();

  doRandomRefinements(adaptivityConfig, *grid, gridGen, alpha);

  sgpp::base::DataVector dataSizeVectorResult(dataset.getNumberInstances());

  eval->prepare();

  eval->mult(alpha, dataSizeVectorResult);

  auto evalCompare = std::shared_ptr<sgpp::base::OperationMultipleEval>(
      sgpp::op_factory::createOperationMultipleEval(*grid, trainingData));

  sgpp::base::DataVector dataSizeVectorResultCompare(dataset.getNumberInstances());

  evalCompare->mult(alpha, dataSizeVectorResultCompare);

  double mse = compareVectors(dataSizeVectorResult, dataSizeVectorResultCompare);

  BOOST_TEST_MESSAGE("fileName: " << fileName << " mse: " << mse);
  return mse;
}

void compareDatasetsTranspose(const std::vector<std::tuple<std::string, double>>& fileNamesError,
                              sgpp::base::GridType gridType, size_t level,
                              sgpp::datadriven::OperationMultipleEvalConfiguration configuration) {
  for (std::tuple<std::string, double> fileNameError : fileNamesError) {
    double mse =
        compareToReferenceTranspose(gridType, std::get<0>(fileNameError), level, configuration);
    BOOST_CHECK(mse < std::get<1>(fileNameError));
    std::cout << "expected error: " << std::get<1>(fileNameError) << ", observed error:" << mse
              << " (transposed)" << std::endl;
  }
}

double compareToReferenceTranspose(
    sgpp::base::GridType gridType, const std::string& fileName, size_t level,
    sgpp::datadriven::OperationMultipleEvalConfiguration configuration) {
  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.numRefinementPoints_ = 80;
  adaptivityConfig.numRefinements_ = 1;
  adaptivityConfig.percent_ = 200.0;
  adaptivityConfig.refinementThreshold_ = 0.0;

  std::string content = uncompressFile(fileName);

  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

  sgpp::base::DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::shared_ptr<sgpp::base::Grid> grid;

  if (gridType == sgpp::base::GridType::Linear) {
    grid = std::shared_ptr<sgpp::base::Grid>(sgpp::base::Grid::createLinearGrid(dim));
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    grid = std::shared_ptr<sgpp::base::Grid>(sgpp::base::Grid::createModLinearGrid(dim));
  }

  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(level);

  sgpp::base::DataVector dataSizeVector(dataset.getNumberInstances());

  // Don't use random data! Random data will change the expected MSE
  for (size_t i = 0; i < dataSizeVector.getSize(); i++) {
    dataSizeVector[i] = static_cast<double>(i + 1);
  }

  auto eval = std::shared_ptr<sgpp::base::OperationMultipleEval>(
      sgpp::op_factory::createOperationMultipleEval(*grid, trainingData, configuration));

  eval->prepare();

  doRandomRefinements(adaptivityConfig, *grid, gridGen);

  sgpp::base::DataVector alphaResult(gridStorage.getSize());

  eval->prepare();

  eval->multTranspose(dataSizeVector, alphaResult);

  auto evalCompare = std::shared_ptr<sgpp::base::OperationMultipleEval>(
      sgpp::op_factory::createOperationMultipleEval(*grid, trainingData));

  sgpp::base::DataVector alphaResultCompare(gridStorage.getSize());

  evalCompare->multTranspose(dataSizeVector, alphaResultCompare);

  double mse = compareVectors(alphaResult, alphaResultCompare);

  BOOST_TEST_MESSAGE("fileName: " << fileName << " mse: " << mse);
  return mse;
}

void compareDatasetsDistributed(const std::vector<std::tuple<std::string, double>>& fileNamesError,
                                sgpp::base::GridType gridType, size_t level,
                                sgpp::datadriven::OperationMultipleEvalConfiguration configuration,
                                std::shared_ptr<BlacsProcessGrid> processGrid) {
  for (std::tuple<std::string, double> fileNameError : fileNamesError) {
    double mse = compareToReferenceDistributed(gridType, std::get<0>(fileNameError), level,
                                               configuration, processGrid);
    BOOST_CHECK(mse < std::get<1>(fileNameError));
    std::cout << "expected error: " << std::get<1>(fileNameError) << ", observed error:" << mse
              << std::endl;
  }
}

double compareToReferenceDistributed(
    sgpp::base::GridType gridType, const std::string& fileName, size_t level,
    sgpp::datadriven::OperationMultipleEvalConfiguration configuration,
    std::shared_ptr<BlacsProcessGrid> processGrid) {
  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.numRefinementPoints_ = 80;
  adaptivityConfig.numRefinements_ = 1;
  adaptivityConfig.percent_ = 200.0;
  adaptivityConfig.refinementThreshold_ = 0.0;

  std::string content = uncompressFile(fileName);

  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

  sgpp::base::DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::shared_ptr<sgpp::base::Grid> grid;

  if (gridType == sgpp::base::GridType::Linear) {
    grid = std::shared_ptr<sgpp::base::Grid>(sgpp::base::Grid::createLinearGrid(dim));
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    grid = std::shared_ptr<sgpp::base::Grid>(sgpp::base::Grid::createModLinearGrid(dim));
  }

  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(level);

  sgpp::base::DataVector alpha(gridStorage.getSize());

  for (size_t i = 0; i < alpha.getSize(); i++) {
    alpha[i] = static_cast<double>(i);
  }

  auto eval = std::shared_ptr<OperationMultipleEvalDistributed>(
      static_cast<OperationMultipleEvalDistributed*>(
          sgpp::op_factory::createOperationMultipleEval(*grid, trainingData, configuration)));

  // no random refinement, as randomness causes problems with multiple processes

  DataVectorDistributed dataSizeVectorResultDistributed(processGrid, dataset.getNumberInstances(),
                                                        32);

  eval->prepare();

  eval->multDistributed(alpha, dataSizeVectorResultDistributed);

  auto evalCompare = std::shared_ptr<sgpp::base::OperationMultipleEval>(
      sgpp::op_factory::createOperationMultipleEval(*grid, trainingData));

  sgpp::base::DataVector dataSizeVectorResultCompare(dataset.getNumberInstances());

  evalCompare->mult(alpha, dataSizeVectorResultCompare);

  auto dataSizeVectorResult = dataSizeVectorResultDistributed.toLocalDataVectorBroadcast();

  double mse = 0;

  if (processGrid->isProcessInGrid()) {
    mse = compareVectors(dataSizeVectorResult, dataSizeVectorResultCompare);
  }

  BOOST_TEST_MESSAGE("fileName: " << fileName << " mse: " << mse);
  return mse;
}

void compareDatasetsTransposeDistributed(
    const std::vector<std::tuple<std::string, double>>& fileNamesError,
    sgpp::base::GridType gridType, size_t level,
    sgpp::datadriven::OperationMultipleEvalConfiguration configuration,
    std::shared_ptr<BlacsProcessGrid> processGrid) {
  for (std::tuple<std::string, double> fileNameError : fileNamesError) {
    double mse = compareToReferenceTransposeDistributed(gridType, std::get<0>(fileNameError), level,
                                                        configuration, processGrid);
    BOOST_CHECK(mse < std::get<1>(fileNameError));
    std::cout << "expected error: " << std::get<1>(fileNameError) << ", observed error:" << mse
              << " (transposed)" << std::endl;
  }
}

double compareToReferenceTransposeDistributed(
    sgpp::base::GridType gridType, const std::string& fileName, size_t level,
    sgpp::datadriven::OperationMultipleEvalConfiguration configuration,
    std::shared_ptr<BlacsProcessGrid> processGrid) {
  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.numRefinementPoints_ = 80;
  adaptivityConfig.numCoarseningPoints_ = 80;
  adaptivityConfig.numRefinements_ = 1;
  adaptivityConfig.percent_ = 200.0;
  adaptivityConfig.refinementThreshold_ = 0.0;

  std::string content = uncompressFile(fileName);

  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

  sgpp::base::DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::shared_ptr<sgpp::base::Grid> grid;

  if (gridType == sgpp::base::GridType::Linear) {
    grid = std::shared_ptr<sgpp::base::Grid>(sgpp::base::Grid::createLinearGrid(dim));
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    grid = std::shared_ptr<sgpp::base::Grid>(sgpp::base::Grid::createModLinearGrid(dim));
  }

  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(level);

  sgpp::base::DataVector dataSizeVector(dataset.getNumberInstances());

  // Don't use random data! Random data will change the expected MSE
  for (size_t i = 0; i < dataSizeVector.getSize(); i++) {
    dataSizeVector[i] = static_cast<double>(i + 1);
  }

  auto eval = std::shared_ptr<OperationMultipleEvalDistributed>(
      static_cast<OperationMultipleEvalDistributed*>(
          sgpp::op_factory::createOperationMultipleEval(*grid, trainingData, configuration)));

  // no random refinement, as randomness causes problems with multiple processes

  DataVectorDistributed alphaResultDistributed(processGrid, gridStorage.getSize(), 32);

  eval->prepare();

  eval->multTransposeDistributed(dataSizeVector, alphaResultDistributed);

  auto evalCompare = std::shared_ptr<sgpp::base::OperationMultipleEval>(
      sgpp::op_factory::createOperationMultipleEval(*grid, trainingData));

  sgpp::base::DataVector alphaResultCompare(gridStorage.getSize());

  evalCompare->multTranspose(dataSizeVector, alphaResultCompare);

  auto alphaResult = alphaResultDistributed.toLocalDataVectorBroadcast();

  double mse = 0;

  if (processGrid->isProcessInGrid()) {
    mse = compareVectors(alphaResult, alphaResultCompare);
  }

  BOOST_TEST_MESSAGE("fileName: " << fileName << " mse: " << mse);
  return mse;
}
#if USE_OCL == 1

std::shared_ptr<OCLOperationConfiguration> getConfigurationDefaultsSingleDevice() {
  // detects the platform
  OCLManagerMultiPlatform manager;
  auto parameters = manager.getConfiguration();

  // filter all devices except for the first
  bool firstPlatform = true;
  bool firstDevice = true;
  std::vector<std::string> platformNames = (*parameters)["PLATFORMS"].keys();
  for (std::string& platformName : platformNames) {
    if (firstPlatform) {
      firstPlatform = false;
    } else {
      (*parameters)["PLATFORMS"][platformName].erase();
      continue;
    }

    json::Node& platformNode = (*parameters)["PLATFORMS"][platformName];

    std::vector<std::string> deviceNames = platformNode["DEVICES"].keys();
    for (std::string& deviceName : deviceNames) {
      if (firstDevice) {
        firstDevice = false;
      } else {
        platformNode["DEVICES"][deviceName].erase();
        continue;
      }

      // make sure there is only a single device of the selected device
      json::Node& deviceNode = platformNode["DEVICES"][deviceName];
      deviceNode.addIDAttr("COUNT", UINT64_C(1));
    }
  }

  return parameters;
}

std::shared_ptr<OCLOperationConfiguration> getConfigurationDefaultsMultiDevice() {
  // detects the platform
  OCLManagerMultiPlatform manager;
  auto parameters = manager.getConfiguration();

  // filter all devices except for the first
  bool firstPlatform = true;
  std::vector<std::string> platformNames = (*parameters)["PLATFORMS"].keys();
  for (std::string& platformName : platformNames) {
    if (firstPlatform) {
      firstPlatform = false;
    } else {
      (*parameters)["PLATFORMS"][platformName].erase();
      continue;
    }
  }

  return parameters;
}

std::shared_ptr<OCLOperationConfiguration> getConfigurationDefaultsMultiPlatform() {
  // detects the platform
  OCLManagerMultiPlatform manager;
  auto parameters = manager.getConfiguration();
  return parameters;
}
#endif
#endif
