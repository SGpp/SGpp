// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/application/LearnerLeastSquaresIdentity.hpp>
#include <sgpp/datadriven/application/MetaLearner.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/globaldef.hpp>

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {

MetaLearner::MetaLearner(sgpp::base::RegularGridConfiguration gridConfig,
                         sgpp::solver::SLESolverConfiguration solverConfig,
                         sgpp::solver::SLESolverConfiguration solverFinalStep,
                         sgpp::base::AdaptivityConfiguration adaptivityConfig, double lambda,
                         bool verbose) {
  this->csvSep = "& ";
  this->gridConfig = gridConfig;
  this->solverConfig = solverConfig;
  this->solverFinalStep = solverFinalStep;
  this->adaptivityConfig = adaptivityConfig;
  this->lambda = lambda;
  this->verbose = verbose;
  this->instances = 0;
}

void MetaLearner::learn(
    sgpp::datadriven::OperationMultipleEvalConfiguration &operationConfiguration,
    std::string &datasetFileName, bool isRegression) {
  std::ifstream t(datasetFileName);
  if (!t.is_open()) {
    throw;
  }
  std::stringstream buffer;
  buffer << t.rdbuf();
  std::string bufferString = buffer.str();
  this->learnString(operationConfiguration, bufferString);
}

void MetaLearner::learnString(
    sgpp::datadriven::OperationMultipleEvalConfiguration &operationConfiguration,
    std::string &content, bool isRegression) {
  Dataset dataset = ARFFTools::readARFFFromString(content);

  this->gridConfig.dim_ = dataset.getDimension();
  this->instances = dataset.getNumberInstances();

  if (verbose) {
    std::cout << "instances: " << this->instances << std::endl;
  }

  base::DataVector &classesVector = dataset.getTargets();
  base::DataMatrix &trainingData = dataset.getData();

  //    bool isRegression = true; // treat everything as if it were a
  //    regression, as classification is not fully supported by Learner
  std::unique_ptr<LearnerLeastSquaresIdentity> myLearner =
      std::make_unique<LearnerLeastSquaresIdentity>(isRegression, this->verbose);
  myLearner->setImplementation(operationConfiguration);
  // TODO(pfandedd): reenabled after performance calculator has been adjusted
  myLearner->setReuseCoefficients(false);

  LearnerTiming timings =
      myLearner->train(trainingData, classesVector, this->gridConfig, this->solverConfig,
                       this->solverFinalStep, this->adaptivityConfig, false, this->lambda);

  this->myTiming = timings;
  this->ExecTimesOnStep = myLearner->getRefinementExecTimes();

  this->myLearner = std::move(myLearner);
}

base::Grid &MetaLearner::getLearnedGrid() {
  if (!this->myLearner.operator bool()) {
    throw base::application_exception("error: cannot get grid if nothing was learned before");
  }
  return this->myLearner->getGrid();
}

base::DataVector &MetaLearner::getLearnedAlpha() {
  if (!this->myLearner.operator bool()) {
    throw base::application_exception("error: cannot get surplusses if nothing was learned before");
  }
  return this->myLearner->getAlpha();
}

void MetaLearner::learnReference(std::string &datasetFileName, bool isRegression) {
  std::ifstream t(datasetFileName);
  if (!t.is_open()) {
    throw;
  }
  std::stringstream buffer;
  buffer << t.rdbuf();
  std::string bufferString = buffer.str();
  this->learnReferenceString(bufferString);
}

void MetaLearner::learnReferenceString(std::string &content, bool isRegression) {
  Dataset dataset = ARFFTools::readARFFFromString(content);
  this->gridConfig.dim_ = dataset.getDimension();
  this->instances = dataset.getNumberInstances();

  if (verbose) {
    std::cout << "instances: " << this->instances << std::endl;
  }

  base::DataVector &classesVector = dataset.getTargets();
  base::DataMatrix &trainingData = dataset.getData();

  std::unique_ptr<LearnerLeastSquaresIdentity> referenceLearner =
      std::make_unique<LearnerLeastSquaresIdentity>(isRegression, this->verbose);
  sgpp::datadriven::OperationMultipleEvalConfiguration operationConfiguration(
      OperationMultipleEvalType::DEFAULT, OperationMultipleEvalSubType::DEFAULT,
      OperationMultipleEvalMPIType::NONE, "STREAMING");
  referenceLearner->setImplementation(operationConfiguration);
  // TODO(pfandedd): reenabled after performance calculator has been adjusted
  referenceLearner->setReuseCoefficients(false);

  LearnerTiming timings =
      referenceLearner->train(trainingData, classesVector, gridConfig, solverConfig,
                              solverFinalStep, adaptivityConfig, false, lambda);
  this->referenceTiming = timings;
  this->ExecTimesOnStepReference = referenceLearner->getRefinementExecTimes();

  this->referenceLearner = std::move(referenceLearner);
}

void MetaLearner::learnAndTest(
    sgpp::datadriven::OperationMultipleEvalConfiguration &operationConfiguration,
    std::string &datasetFileName, std::string &testFileName, bool isRegression) {
  std::ifstream dataFile(datasetFileName);
  std::stringstream bufferData;
  bufferData << dataFile.rdbuf();
  std::ifstream testFile(datasetFileName);
  std::stringstream bufferTest;
  bufferTest << testFile.rdbuf();
  std::string bufferDataString = bufferData.str();
  std::string bufferTestString = bufferTest.str();
  this->learnAndTestString(operationConfiguration, bufferDataString, bufferTestString,
                           isRegression);
}

// learn and test against test dataset and measure hits/mse
void MetaLearner::learnAndTestString(
    sgpp::datadriven::OperationMultipleEvalConfiguration &operationConfiguration,
    std::string &dataContent, std::string &testContent, bool isRegression) {
  // always to this first
  this->learnString(operationConfiguration, dataContent);

  Dataset testDataset = ARFFTools::readARFFFromString(testContent);
  size_t testDim = testDataset.getDimension();
  size_t testInstances = testDataset.getNumberInstances();

  base::DataVector &testClassesVector = testDataset.getTargets();
  base::DataMatrix &testTrainingData = testDataset.getData();

  if (verbose && testDim != this->gridConfig.dim_) {
    std::cout << "dim of test dataset and training dataset doesn't match" << std::endl;
  }

  if (verbose) {
    std::cout << "computing classes of test dataset" << std::endl;
  }

  base::DataVector computedClasses(testTrainingData.getNrows());
  this->myLearner->predict(testTrainingData, computedClasses);

  if (verbose) {
    std::cout << "classes computed" << std::endl;
  }

  if (isRegression) {
    double mse = 0.0;

    for (size_t i = 0; i < computedClasses.getSize(); i++) {
      double diff = testClassesVector.get(i) - computedClasses.get(i);
      mse += diff * diff;
    }

    mse = mse / static_cast<double>(testInstances);

    if (verbose) {
      std::cout << "mse: " << mse << std::endl;
    }
  } else {
    int hits = 0;

    for (size_t i = 0; i < computedClasses.getSize(); i++) {
      int classToken = 0;

      if (computedClasses.get(i) > 0) {
        classToken = 1;
      } else {
        classToken = -1;
      }

      if (classToken == testClassesVector.get(i)) {
        hits += 1;
      }
    }

    if (verbose) {
      std::cout << "hits (%): " << (static_cast<double>(hits) / static_cast<double>(testInstances))
                << std::endl;
    }
  }
}

// learn and test against the streaming implemenation
double MetaLearner::learnAndCompare(
    sgpp::datadriven::OperationMultipleEvalConfiguration &operationConfiguration,
    std::string &datasetFileName, size_t gridGranularity) {
  std::ifstream t(datasetFileName);
  if (!t.is_open()) {
    throw;
  }
  std::stringstream buffer;
  buffer << t.rdbuf();
  std::string bufferString = buffer.str();
  return this->learnAndCompareString(operationConfiguration, bufferString, gridGranularity);
}

// learn and test against the streaming implemenation
double MetaLearner::learnAndCompareString(
    sgpp::datadriven::OperationMultipleEvalConfiguration &operationConfiguration,
    std::string &content, size_t gridGranularity) {
  // always do this first
  this->learnString(operationConfiguration, content);
  this->learnReferenceString(content);

  base::DataMatrix testTrainingData(0, this->gridConfig.dim_);

  double increment = 1.0 / static_cast<double>(gridGranularity);

  size_t index = 0;
  base::DataVector testPoint(this->gridConfig.dim_);

  for (size_t i = 0; i < this->gridConfig.dim_; i++) {
    testPoint[i] = 0.0 + increment;
  }

  testTrainingData.appendRow(testPoint);

  unsigned int testInstanceCounter = 0;

  while (index < this->gridConfig.dim_) {
    // make 1.0 impossible
    double prior = testPoint[index] + increment;

    if (prior < 1.0) {
      testPoint[index] += increment;

      for (size_t i = 0; i < index; i++) {
        testPoint[i] = 0.0 + increment;  // skip border
      }

      testTrainingData.appendRow(testPoint);
      testInstanceCounter += 1;

      if (verbose && testInstanceCounter % 1000000 == 0) {
        std::cout << "testInstanceCounter (still generating): " << testInstanceCounter << std::endl;
      }

      index = 0;
    } else {
      index += 1;
    }
  }

  if (verbose) {
    std::cout << "testInstanceCounter: " << (testInstanceCounter + 1) << std::endl;
    std::cout << "predicting..." << std::endl;
  }

  base::DataVector computedClasses(testTrainingData.getNrows());
  this->myLearner->predict(testTrainingData, computedClasses);

  if (verbose) {
    std::cout << "predicting... (reference)" << std::endl;
  }

  base::DataVector referenceClasses(testTrainingData.getNrows());
  this->referenceLearner->predict(testTrainingData, referenceClasses);

  double squareSum = 0.0;

  // sgpp::base::DataVector *myAlpha = this->myLearner->alpha_;
  // for (size_t i = 0; i < myAlpha->getSize();i++) {
  //    std::cout << "alpha[ " << i << "]=" << (*myAlpha)[i] << ", ";
  // }
  // std::cout << std::endl;

  double maxDiff = -1.0;
  double firstValue = -1.0;
  double secondValue = -1.0;
  for (size_t i = 0; i < computedClasses.getSize(); i++) {
    double diff = fabs(computedClasses.get(i) - referenceClasses.get(i));

    if (diff > maxDiff) {
      maxDiff = diff;
      firstValue = computedClasses.get(i);
      secondValue = referenceClasses.get(i);
      //      std::cout << "first: " << firstValue << " second: " << secondValue << std::endl;
    }

    diff *= diff;
    squareSum += diff;
  }

  squareSum = sqrt(squareSum);

  if (verbose) {
    std::cout << "sqrt: " << squareSum << std::endl;
    std::cout << "maxDiff: " << maxDiff << " value: " << firstValue << " second: " << secondValue
              << std::endl;
  }

  return squareSum;
}

LearnerTiming MetaLearner::getLearnerTiming() { return this->myTiming; }

LearnerTiming MetaLearner::getLearnerReferenceTiming() { return this->referenceTiming; }
}  // namespace datadriven
}  // namespace sgpp
