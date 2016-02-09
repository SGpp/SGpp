// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/MetaLearner.hpp>

#include <sgpp/datadriven/application/LearnerLeastSquaresIdentity.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/tools/DatasetTools.hpp>

#include <sgpp/globaldef.hpp>

#include <random>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <vector>

namespace SGPP {
namespace datadriven {

MetaLearner::MetaLearner(SGPP::base::RegularGridConfiguration gridConfig,
                         SGPP::solver::SLESolverConfiguration solverConfig,
                         SGPP::solver::SLESolverConfiguration solverFinalStep,
                         SGPP::base::AdpativityConfiguration adaptivityConfiguration,
                         bool verbose) {
  this->csvSep = "& ";
  this->gridConfig = gridConfig;
  this->solverConfig = solverConfig;
  this->solverFinalStep = solverFinalStep;
  this->adaptivityConfiguration = adaptivityConfiguration;
  this->verbose = verbose;
  this->instances = 0;
}

void MetaLearner::learn(SGPP::datadriven::OperationMultipleEvalConfiguration&
                        operationConfiguration,
                        std::string& datasetFileName, float_t lambda, bool isRegression) {
  std::ifstream t(datasetFileName);

  if (!t.is_open()) {
    throw;
  }

  std::stringstream buffer;
  buffer << t.rdbuf();
  std::string bufferString = buffer.str();
  this->learnString(operationConfiguration, bufferString, lambda, isRegression);
}

void MetaLearner::learnString(
    SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
    std::string& content, float_t lambda, bool isRegression) {
  Dataset dataset = ARFFTools::readARFFFromString(content);

  this->gridConfig.dim_ = dataset.getDimension();
  this->instances = dataset.getNumberInstances();

  if (verbose) {
    std::cout << "instances: " << this->instances << std::endl;
  }

  base::DataVector& classesVector = dataset.getTargets();
  base::DataMatrix& trainingData = dataset.getData();

  //    bool isRegression = true; // treat everything as if it were a regression, as classification
  // is not fully supported by Learner
  auto myLearner = std::make_shared<LearnerLeastSquaresIdentity>(isRegression,
                   this->verbose);
  myLearner->setImplementation(operationConfiguration);

  LearnerTiming timings = myLearner->train(trainingData, classesVector,
                          this->gridConfig, this->solverConfig,
                          this->solverFinalStep, this->adaptivityConfiguration, false, lambda);

  this->myTiming = timings;
  this->ExecTimesOnStep = myLearner->getRefinementExecTimes();

  this->myLearner = myLearner;
}

std::shared_ptr<base::Grid> MetaLearner::getLearnedGrid() {
  if (this->myLearner == nullptr) {
    throw;
  }

  return this->myLearner->getGridCopy();
}

void MetaLearner::learnReference(std::string& datasetFileName, float_t lambda,
                                 bool isRegression) {
  std::ifstream t(datasetFileName);

  if (!t.is_open()) {
    throw;
  }

  std::stringstream buffer;
  buffer << t.rdbuf();
  std::string bufferString = buffer.str();
  this->learnReferenceString(bufferString, lambda, isRegression);
}

void MetaLearner::learnReferenceString(std::string& content, float_t lambda,
                                       bool isRegression) {
  Dataset dataset = ARFFTools::readARFFFromString(content);
  this->gridConfig.dim_ = dataset.getDimension();
  this->instances = dataset.getNumberInstances();

  if (verbose) {
    std::cout << "instances: " << this->instances << std::endl;
  }

  base::DataVector& classesVector = dataset.getTargets();
  base::DataMatrix& trainingData = dataset.getData();

  // treat everything as if it were a regression, as classification is not fully
  // supported by Learner
  //    bool isRegression = true;
  auto referenceLearner = std::make_shared<LearnerLeastSquaresIdentity>
                          (isRegression, this->verbose);
  SGPP::datadriven::OperationMultipleEvalConfiguration operationConfiguration(
    OperationMultipleEvalType::DEFAULT,
    OperationMultipleEvalSubType::DEFAULT, "STREAMING");
  referenceLearner->setImplementation(operationConfiguration);

  LearnerTiming timings = referenceLearner->train(trainingData, classesVector,
                          gridConfig, solverConfig,
                          solverFinalStep, adaptivityConfiguration, false, lambda);
  this->referenceTiming = timings;
  this->ExecTimesOnStepReference = referenceLearner->getRefinementExecTimes();

  // referenceLearner->dumpFunction("referenceGridFunction", 50);
  this->referenceLearner = referenceLearner;
}

void MetaLearner::learnAndTest(
  SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
  std::string& datasetFileName, std::string& testFileName, bool isRegression) {
  std::ifstream dataFile(datasetFileName);
  std::stringstream bufferData;
  bufferData << dataFile.rdbuf();
  std::ifstream testFile(datasetFileName);
  std::stringstream bufferTest;
  bufferTest << testFile.rdbuf();
  std::string bufferDataString = bufferData.str();
  std::string bufferTestString = bufferTest.str();
  this->learnAndTestString(operationConfiguration, bufferDataString,
                           bufferTestString, isRegression);
}

// learn and test against test dataset and measure hits/mse
void MetaLearner::learnAndTestString(
  SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
  std::string& dataContent, std::string& testContent, float_t lambda,
  bool isRegression) {
  // always to this first
  this->learnString(operationConfiguration, dataContent, lambda);

  Dataset testDataset = ARFFTools::readARFFFromString(testContent);
  size_t testDim = testDataset.getDimension();
  size_t testInstances = testDataset.getNumberInstances();

  base::DataVector& testClassesVector = testDataset.getTargets();
  base::DataMatrix& testTrainingData = testDataset.getData();

  if (verbose && testDim != this->gridConfig.dim_) {
    std::cout << "dim of test dataset and training dataset doesn't match" <<
              std::endl;
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
    float_t mse = 0.0;

    for (size_t i = 0; i < computedClasses.getSize(); i++) {
      float_t diff = testClassesVector.get(i) - computedClasses.get(i);
      mse += diff * diff;
    }

    mse = mse / (float_t) testInstances;

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
      std::cout << "hits (%): " << ((float_t) hits / (float_t) testInstances) <<
                std::endl;
    }
  }
}

// learn and test against the streaming implemenation
float_t MetaLearner::learnAndCompare(
  SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
  std::string& datasetFileName, float_t lambda, size_t gridGranularity) {
  std::ifstream t(datasetFileName);

  if (!t.is_open()) {
    throw;
  }

  std::stringstream buffer;
  buffer << t.rdbuf();
  std::string bufferString = buffer.str();
  return this->learnAndCompareString(operationConfiguration, bufferString, lambda,
                                     gridGranularity);
}

// learn and test against the streaming implemenation
float_t MetaLearner::learnAndCompareString(
  SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
  std::string& content, float_t lambda, size_t gridGranularity) {
  // always do this first
  this->learnString(operationConfiguration, content, lambda);
  this->learnReferenceString(content, lambda);

  base::DataMatrix testTrainingData(0, this->gridConfig.dim_);

  float_t increment = 1.0 / static_cast<float_t>(gridGranularity);

  size_t index = 0;
  base::DataVector testPoint(this->gridConfig.dim_);

  for (size_t i = 0; i < this->gridConfig.dim_; i++) {
    testPoint[i] = 0.0 + increment;
  }

  testTrainingData.appendRow(testPoint);

  unsigned int testInstanceCounter = 0;

  while (index < this->gridConfig.dim_) {
    // make 1.0 impossible
    float_t prior = testPoint[index] + increment;

    if (prior < 1.0) {
      testPoint[index] += increment;

      for (size_t i = 0; i < index; i++) {
        testPoint[i] = 0.0 + increment;  // skip border
      }

      testTrainingData.appendRow(testPoint);
      testInstanceCounter += 1;

      if (verbose && testInstanceCounter % 1000000 == 0) {
        std::cout << "testInstanceCounter (still generating): " << testInstanceCounter
                  << std::endl;
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

  float_t squareSum = 0.0;

  double maxDiff = -1.0;
  double firstValue = -1.0;
  double secondValue = -1.0;

  for (size_t i = 0; i < computedClasses.getSize(); i++) {
    float_t diff = fabs(computedClasses.get(i) - referenceClasses.get(i));

    if (diff > maxDiff) {
      maxDiff = diff;
      firstValue = computedClasses.get(i);
      secondValue = referenceClasses.get(i);
    }

    diff *= diff;
    squareSum += diff;
  }

  squareSum = sqrt(squareSum);

  if (verbose) {
    std::cout << "sqrt: " << squareSum << std::endl;
    std::cout << "maxDiff: " << maxDiff << " value: " << firstValue << " second: "
              << secondValue << std::endl;
  }

  return squareSum;
}

LearnerTiming MetaLearner::getLearnerTiming() {
  return this->myTiming;
}

LearnerTiming MetaLearner::getLearnerReferenceTiming() {
  return this->referenceTiming;
}

void MetaLearner::optimizeLambdaLog(SGPP::base::DataMatrix& dataset,
                                    SGPP::base::DataVector& values, size_t kFold,
                                    size_t maxLevel, std::shared_ptr<base::Grid>& bestGrid,
                                    std::shared_ptr<base::DataVector>& bestAlpha,
                                    float_t& lambdaOpt, datadriven::LearnerTiming& timing) {
  lambdaOpt = this->optimizeLambdaLog(dataset, values, kFold, maxLevel);
  this->train(dataset, values, lambdaOpt, bestGrid, bestAlpha, timing);
}

float_t MetaLearner::optimizeLambdaLog(SGPP::base::DataMatrix& dataset,
                                       SGPP::base::DataVector& values, size_t kFold,
                                       size_t maxLevel) {
  std::vector<base::DataMatrix> trainingSets;
  std::vector<base::DataVector> trainingSetsValues;
  std::vector<base::DataMatrix> testSets;
  std::vector<base::DataVector> testSetsValues;

  datadriven::DatasetTools::splitset(dataset, values, kFold, trainingSets,
                                     trainingSetsValues, testSets, testSetsValues);

  // initial values are pure dummy values
  float_t bestLambda = 0.0;
  float_t bestMSE = 0.0;

  this->optimizeLambdaLog_(dataset, values, kFold, maxLevel, trainingSets,
                           trainingSetsValues, testSets,
                           testSetsValues, 0, 1.0, bestLambda, bestMSE);
  bestLambda = pow(10.0, -bestLambda);

  if (verbose) {
    std::cout << "# -> bestLambda = " << bestLambda << std::endl;
    std::cout << "# -> bestMSE= " << bestMSE << std::endl;
  }

  return bestLambda;
}

void MetaLearner::optimizeLambdaLog_(SGPP::base::DataMatrix& dataset,
                                     SGPP::base::DataVector& datasetValues,
                                     size_t kFold, size_t maxLevel,
                                     std::vector<base::DataMatrix>& trainingSets,
                                     std::vector<base::DataVector>& trainingSetsValues,
                                     std::vector<base::DataMatrix>& testSets,
                                     std::vector<base::DataVector>& testSetsValues,
                                     size_t curLevel,
                                     float_t lambdaLogStepSize,
                                     float_t& bestLogLambda, float_t& bestMSE) {
  if (verbose) {
    std::cout << "entering level=" << curLevel << " with lambda=" << pow(10.0,
              -bestLogLambda) << std::endl;
  }

  std::vector<float_t> logLambdaValues;

  if (curLevel == 0 && lambdaLogStepSize == 1.0) {
    for (size_t i = 3; i <= 10; i++) {
      logLambdaValues.push_back(static_cast<float_t>(i));
    }
  } else {
    logLambdaValues.push_back(bestLogLambda - lambdaLogStepSize);

    if (curLevel == 0) {
      logLambdaValues.push_back(bestLogLambda);
    }

    logLambdaValues.push_back(bestLogLambda + lambdaLogStepSize);
  }

  bool firstValue = true;

  for (float_t curLogLambda : logLambdaValues) {
    float_t curLambda = pow(10, -curLogLambda);
    std::cout << "curLambda: " << curLambda << std::endl;

    // cross-validation
    float_t curMeanMSE = 0.0;

    for (size_t j = 0; j < kFold; j++) {
      std::shared_ptr<base::Grid> grid;
      std::shared_ptr<base::DataVector> alpha;
      datadriven::LearnerTiming timing;

      // compute density
      train(trainingSets[j], trainingSetsValues[j], curLambda, grid, alpha, timing);

      float_t mse = this->calculateMSE(*grid, *alpha, testSets[j], testSetsValues[j]);
      curMeanMSE += mse;
    }

    curMeanMSE /= static_cast<float_t>(kFold);

    if ((curLevel == 0 && firstValue) || curMeanMSE < bestMSE) {
      bestMSE = curMeanMSE;
      bestLogLambda = curLogLambda;
      firstValue = false;

      if (verbose) {
        std::cout << "new best lambda!" << std::endl;
      }
    } else {
      if (curLevel == 0) {
        break;
      }
    }

    if (verbose) {
      std::cout << "# lambda: " << curLambda << " curMeanMSE: " << curMeanMSE <<
                " bestLambda: "
                << pow(10.0, -bestLogLambda) << " bestMSE: " << bestMSE <<
                " lambdaLogStepSize: "
                << lambdaLogStepSize << std::endl;
    }
  }

  if (curLevel < maxLevel) {
    this->optimizeLambdaLog_(dataset, datasetValues, kFold, maxLevel, trainingSets,
                             trainingSetsValues, testSets,
                             testSetsValues, curLevel + 1, lambdaLogStepSize / 2.0,
                             bestLogLambda, bestMSE);
  }
}

void MetaLearner::train(base::DataMatrix& train, base::DataVector& trainValues,
                        float_t lambda,
                        std::shared_ptr<base::Grid>& grid, std::shared_ptr<base::DataVector>& alpha,
                        datadriven::LearnerTiming& timing) {
  if (verbose) {
    std::cout << "training..." << std::endl;
  }

  bool isRegression = true;

  auto myLearner = std::make_shared<LearnerLeastSquaresIdentity>(isRegression,
                   this->verbose);

  SGPP::datadriven::OperationMultipleEvalConfiguration config(
    OperationMultipleEvalType::DEFAULT,
    OperationMultipleEvalSubType::DEFAULT);
  myLearner->setImplementation(config);

  timing = myLearner->train(train, trainValues, this->gridConfig,
                            this->solverConfig, this->solverFinalStep,
                            this->adaptivityConfiguration, false, lambda);

  this->myLearner = myLearner;

  grid = std::shared_ptr<base::Grid>(myLearner->getGridCopy());
  alpha = std::shared_ptr<base::DataVector>(myLearner->getAlphaCopy());

  if (verbose) {
    std::cout << "training finished" << std::endl;
  }
}

float_t MetaLearner::calculateMSE(base::Grid& grid, base::DataVector& alpha,
                                  base::DataMatrix& testSubset,
                                  base::DataVector& valuesTestSubset, bool verbose) {
  size_t dim = testSubset.getNcols();
  float_t mse = 0.0;

  base::OperationEval* opEval = SGPP::op_factory::createOperationEval(grid);

  for (size_t i = 0; i < testSubset.getNrows(); i++) {
    base::DataVector point(dim);
    testSubset.getRow(i, point);
    float_t approximation = opEval->eval(alpha, point);
    mse += (approximation - valuesTestSubset[i]) * (approximation -
           valuesTestSubset[i]);

    if (verbose && i < 100) {
      std::cout << "mine: " << approximation << " reference: " << valuesTestSubset[i]
                << std::endl;
    }
  }

  return mse;
}

}  // namespace datadriven
}  // namespace SGPP

