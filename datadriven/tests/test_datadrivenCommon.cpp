// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <zlib.h>
#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <random>
#include <fstream>

#include "test_datadrivenCommon.hpp"
#include "sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"

using SGPP::base::DataMatrix;
using SGPP::base::DataVector;
#if USE_OCL == 1
#include "sgpp/base/opencl/OCLManagerMultiPlatform.hpp"
using SGPP::base::OCLManagerMultiPlatform;
using SGPP::base::OCLOperationConfiguration;
#endif

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

DataMatrix* readReferenceMatrix(SGPP::base::GridStorage& storage, std::string fileName) {
  std::string content = uncompressFile(fileName);

  std::stringstream contentStream;
  contentStream << content;
  std::string line;

  DataMatrix* m = new DataMatrix(0, storage.size());

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
    float_t floatValue;

    for (size_t i = 0; i < storage.size(); i++) {
      curFind = line.find_first_of(" \t", curPos);
      curValue = line.substr(curPos, curFind - curPos);

      floatValue = boost::lexical_cast<float_t>(curValue);
      m->set(currentRow, i, floatValue);
      curPos = curFind + 1;
    }

    currentRow += 1;
  }

  return m;
}

void doRandomRefinements(SGPP::base::AdpativityConfiguration& adaptConfig, SGPP::base::Grid& grid,
                         SGPP::base::GridGenerator& gridGen, SGPP::base::DataVector& alpha) {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  for (size_t i = 0; i < adaptConfig.numRefinements_; i++) {
    SGPP::base::SurplusRefinementFunctor myRefineFunc(
        alpha, adaptConfig.noPoints_, adaptConfig.threshold_);
    gridGen.refine(myRefineFunc);
    size_t oldSize = alpha.getSize();
    alpha.resize(grid.getSize());

    for (size_t j = oldSize; j < alpha.getSize(); j++) {
      alpha[j] = dist(mt);
    }
  }
}

void doRandomRefinements(SGPP::base::AdpativityConfiguration& adaptConfig, SGPP::base::Grid& grid,
                         SGPP::base::GridGenerator& gridGen) {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  SGPP::base::DataVector alphaRefine(grid.getSize());

  for (size_t i = 0; i < alphaRefine.getSize(); i++) {
    alphaRefine[i] = dist(mt);
  }

  for (size_t i = 0; i < adaptConfig.numRefinements_; i++) {
    SGPP::base::SurplusRefinementFunctor myRefineFunc(
        alphaRefine, adaptConfig.noPoints_, adaptConfig.threshold_);
    gridGen.refine(myRefineFunc);
    size_t oldSize = alphaRefine.getSize();
    alphaRefine.resize(grid.getSize());

    for (size_t j = oldSize; j < alphaRefine.getSize(); j++) {
      alphaRefine[j] = dist(mt);
    }
  }
}

double compareVectors(SGPP::base::DataVector& results, SGPP::base::DataVector& resultsCompare) {
  double mse = 0.0;

  bool anyDifferentValue = false;
  double largestDifference = 0.0;
  double value = 0.0;
  double valueReference = 0.0;

  for (size_t i = 0; i < resultsCompare.getSize(); i++) {
    double diff = (results[i] - resultsCompare[i]) * (results[i] - resultsCompare[i]);

    if (diff > largestDifference) {
      anyDifferentValue = true;
      largestDifference = diff;
      value = results[i];
      valueReference = resultsCompare[i];
    }

    //        BOOST_TEST_MESSAGE("i: " << i << " mine: " << results[i] << " ref: " <<
    //        resultsCompare[i]);
    mse += (results[i] - resultsCompare[i]) * (results[i] - resultsCompare[i]);
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

double compareToReference(SGPP::base::GridType gridType, std::string fileName, size_t level,
                          SGPP::datadriven::OperationMultipleEvalConfiguration configuration,
                          size_t numRefinements) {
  SGPP::base::AdpativityConfiguration adaptConfig;
  adaptConfig.maxLevelType_ = false;
  adaptConfig.noPoints_ = 80;
  adaptConfig.numRefinements_ = numRefinements;
  adaptConfig.percent_ = 200.0;
  adaptConfig.threshold_ = 0.0;

  std::string content = uncompressFile(fileName);

  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

  SGPP::base::DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::shared_ptr<SGPP::base::Grid> grid;

  if (gridType == SGPP::base::GridType::Linear) {
    grid = std::shared_ptr<SGPP::base::Grid>(SGPP::base::Grid::createLinearGrid(dim));
  } else if (gridType == SGPP::base::GridType::ModLinear) {
    grid = std::shared_ptr<SGPP::base::Grid>(SGPP::base::Grid::createModLinearGrid(dim));
  }

  SGPP::base::GridStorage& gridStorage = grid->getStorage();

  auto gridGen = std::shared_ptr<SGPP::base::GridGenerator>(grid->createGridGenerator());
  gridGen->regular(level);

  SGPP::base::DataVector alpha(gridStorage.size());

  for (size_t i = 0; i < alpha.getSize(); i++) {
    alpha[i] = static_cast<double>(i);
  }

  auto eval = std::shared_ptr<SGPP::base::OperationMultipleEval>(
      SGPP::op_factory::createOperationMultipleEval(*grid, trainingData, configuration));

  eval->prepare();

  doRandomRefinements(adaptConfig, *grid, *gridGen, alpha);

  SGPP::base::DataVector dataSizeVectorResult(dataset.getNumberInstances());
  dataSizeVectorResult.setAll(0);

  eval->prepare();

  eval->mult(alpha, dataSizeVectorResult);

  auto evalCompare = std::shared_ptr<SGPP::base::OperationMultipleEval>(
      SGPP::op_factory::createOperationMultipleEval(*grid, trainingData));

  SGPP::base::DataVector dataSizeVectorResultCompare(dataset.getNumberInstances());
  dataSizeVectorResultCompare.setAll(0.0);

  evalCompare->mult(alpha, dataSizeVectorResultCompare);

  double mse = compareVectors(dataSizeVectorResult, dataSizeVectorResultCompare);

  BOOST_TEST_MESSAGE("fileName: " << fileName << " mse: " << mse);
  return mse;
}

double compareToReferenceTranspose(
    SGPP::base::GridType gridType, std::string fileName, size_t level,
    SGPP::datadriven::OperationMultipleEvalConfiguration configuration) {
  SGPP::base::AdpativityConfiguration adaptConfig;
  adaptConfig.maxLevelType_ = false;
  adaptConfig.noPoints_ = 80;
  adaptConfig.numRefinements_ = 1;
  adaptConfig.percent_ = 200.0;
  adaptConfig.threshold_ = 0.0;

  std::string content = uncompressFile(fileName);

  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

  SGPP::base::DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::shared_ptr<SGPP::base::Grid> grid;

  if (gridType == SGPP::base::GridType::Linear) {
    grid = std::shared_ptr<SGPP::base::Grid>(SGPP::base::Grid::createLinearGrid(dim));
  } else if (gridType == SGPP::base::GridType::ModLinear) {
    grid = std::shared_ptr<SGPP::base::Grid>(SGPP::base::Grid::createModLinearGrid(dim));
  }

  SGPP::base::GridStorage& gridStorage = grid->getStorage();

  auto gridGen = std::shared_ptr<SGPP::base::GridGenerator>(grid->createGridGenerator());
  gridGen->regular(level);

  SGPP::base::DataVector dataSizeVector(dataset.getNumberInstances());

  // Don't use random data! Random data will change the expected MSE
  for (size_t i = 0; i < dataSizeVector.getSize(); i++) {
    dataSizeVector[i] = static_cast<double>(i + 1);
  }

  auto eval = std::shared_ptr<SGPP::base::OperationMultipleEval>(
      SGPP::op_factory::createOperationMultipleEval(*grid, trainingData, configuration));

  eval->prepare();

  doRandomRefinements(adaptConfig, *grid, *gridGen);

  SGPP::base::DataVector alphaResult(gridStorage.size());

  eval->prepare();

  eval->multTranspose(dataSizeVector, alphaResult);

  auto evalCompare = std::shared_ptr<SGPP::base::OperationMultipleEval>(
      SGPP::op_factory::createOperationMultipleEval(*grid, trainingData));

  SGPP::base::DataVector alphaResultCompare(gridStorage.size());
  alphaResultCompare.setAll(0.0);

  evalCompare->multTranspose(dataSizeVector, alphaResultCompare);

  double mse = compareVectors(alphaResult, alphaResultCompare);

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
  for (std::string& platformName : (*parameters)["PLATFORMS"].keys()) {
    if (firstPlatform) {
      firstPlatform = false;
    } else {
      (*parameters)["PLATFORMS"][platformName].erase();
      continue;
    }

    json::Node& platformNode = (*parameters)["PLATFORMS"][platformName];

    for (std::string& deviceName : platformNode["DEVICES"].keys()) {
      if (firstDevice) {
        firstDevice = false;
      } else {
        platformNode["DEVICES"][deviceName].erase();
        continue;
      }

      // make sure there is only a single device of the selected device
      json::Node& deviceNode = platformNode["DEVICES"][deviceName];
      deviceNode.addIDAttr("COUNT", 1ul);
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
  for (std::string& platformName : (*parameters)["PLATFORMS"].keys()) {
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
