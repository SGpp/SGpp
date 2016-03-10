// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "MetaLearner.hpp"

#include <random>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <utility>
#include <string>
#include <vector>

#include "sgpp/globaldef.hpp"
#include "LearnerLeastSquaresIdentity.hpp"
#include "sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/datadriven/tools/Dataset.hpp"

namespace sgpp {
namespace datadriven {

MetaLearner::MetaLearner(sgpp::base::RegularGridConfiguration gridConfig,
                         sgpp::solver::SLESolverConfiguration solverConfig,
                         sgpp::solver::SLESolverConfiguration solverFinalStep,
                         sgpp::base::AdpativityConfiguration adaptivityConfiguration,
                         double lambda, bool verbose) {
  this->csvSep = "& ";
  this->gridConfig = gridConfig;
  this->solverConfig = solverConfig;
  this->solverFinalStep = solverFinalStep;
  this->adaptivityConfiguration = adaptivityConfiguration;
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
  LearnerLeastSquaresIdentity *myLearner =
      new LearnerLeastSquaresIdentity(isRegression, this->verbose);
  myLearner->setImplementation(operationConfiguration);

  LearnerTiming timings =
      myLearner->train(trainingData, classesVector, this->gridConfig, this->solverConfig,
                       this->solverFinalStep, this->adaptivityConfiguration, false, this->lambda);

  this->myTiming = timings;
  this->ExecTimesOnStep = myLearner->getRefinementExecTimes();

  if (this->myLearner != nullptr) {
    delete this->myLearner;
  }
  this->myLearner = myLearner;
}

base::Grid &MetaLearner::getLearnedGrid() {
  if (this->myLearner == nullptr) {
    throw;
  }
  return this->myLearner->getGrid();
}

base::DataVector &MetaLearner::getLearnedAlpha() {
  if (this->myLearner == nullptr) {
    throw;
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

  //    // treat everything as if it were a regression, as classification is not
  //    fully supported by Learner
  //    bool isRegression = true;
  LearnerLeastSquaresIdentity *referenceLearner =
      new LearnerLeastSquaresIdentity(isRegression, this->verbose);
  sgpp::datadriven::OperationMultipleEvalConfiguration operationConfiguration(
      OperationMultipleEvalType::DEFAULT, OperationMultipleEvalSubType::DEFAULT, "STREAMING");
  referenceLearner->setImplementation(operationConfiguration);

  LearnerTiming timings =
      referenceLearner->train(trainingData, classesVector, gridConfig, solverConfig,
                              solverFinalStep, adaptivityConfiguration, false, lambda);
  this->referenceTiming = timings;
  this->ExecTimesOnStepReference = referenceLearner->getRefinementExecTimes();

  // referenceLearner->dumpFunction("referenceGridFunction", 50);
  if (this->referenceLearner != nullptr) {
    delete this->referenceLearner;
  }

  this->referenceLearner = referenceLearner;
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
      std::cout << "hits (%): " << (static_cast<double>(hits) /
          static_cast<double>(testInstances)) << std::endl;
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

void MetaLearner::writeRefinementResults(
    std::string fileName, std::string fileHeader,
    std::vector<std::pair<std::string, std::vector<std::pair<size_t, double> > > > datasetDetails,
    std::vector<std::pair<std::string, std::vector<std::pair<size_t, double> > > >
        datasetDetailsReference,
    bool referenceComparison) {
  std::ofstream myFile;
  myFile.open(fileName, std::ios::out);
  myFile << "# " << fileHeader << std::endl;
  myFile << "refinement";

  for (auto datasetDetail : datasetDetails) {
    myFile << csvSep << datasetDetail.first << "/mine";

    if (referenceComparison) {
      myFile << csvSep << datasetDetail.first << "/reference";
      myFile << csvSep << datasetDetail.first << "/speedup";
    }
  }

  myFile << std::endl;
  size_t refinements = datasetDetails[0].second.size();

  for (size_t i = 0; i < refinements; i++) {
    bool first = true;

    for (size_t datasetIndex = 0; datasetIndex < datasetDetails.size(); datasetIndex++) {
      std::pair<size_t, double> stepTuple = datasetDetails[datasetIndex].second[i];
      double executionTime = stepTuple.second;

      if (first) {
        // refinement steps
        myFile << stepTuple.first;
        first = false;
      }

      // execution time at refinement step
      myFile << csvSep << executionTime;

      if (referenceComparison) {
        std::pair<std::string, std::vector<std::pair<size_t, double> > > datasetTuple =
            datasetDetailsReference[datasetIndex];
        std::pair<size_t, double> stepTupleReference = datasetTuple.second[i];
        double executionTimeReference = stepTupleReference.second;
        myFile << csvSep << executionTimeReference;
        // speedup
        myFile << csvSep << executionTimeReference / executionTime;
      }
    }

    myFile << std::endl;
  }

  myFile.close();
}

void MetaLearner::refinementAndOverallPerformance(
    std::vector<sgpp::datadriven::OperationMultipleEvalConfiguration *> operationConfigurations,
    std::vector<std::string> datasets, std::vector<std::string> experimentHeaders,
    std::string metaInformation, std::string experimentName, bool referenceComparison) {
  for (sgpp::datadriven::OperationMultipleEvalConfiguration *operationConfiguration :
       operationConfigurations) {
    std::vector<std::pair<std::string, std::vector<std::pair<size_t, double> > > >
        refinementDetails;
    std::vector<std::pair<std::string, std::vector<std::pair<size_t, double> > > >
        refinementDetailsReference;

    std::string kernelMetaInformation =
        metaInformation + " kernel: " + operationConfiguration->getName();

    std::ofstream myFile;
    std::string experimentFile =
        "results/data/overall_" + experimentName + "_" + operationConfiguration->getName() + ".csv";
    remove(experimentFile.c_str());
    std::string experimentDetailsFile = "results/data/refinement_" + experimentName + "_" +
                                        operationConfiguration->getName() + ".csv";
    myFile.open(experimentFile, std::ios::out | std::ios::app);
    myFile << "# " << kernelMetaInformation << std::endl;
    myFile << "experiment" << csvSep << "duration" << csvSep << "durationLastIt";

    if (referenceComparison) {
      myFile << csvSep << "durationRef";
      myFile << csvSep << "durationLastItRef";
      myFile << csvSep << "speedup";
      myFile << csvSep << "speedupLastIt";
    }

    myFile << std::endl;

    for (size_t i = 0; i < datasets.size(); i++) {
      this->learn(*operationConfiguration, datasets[i]);
      refinementDetails.push_back(make_pair(datasets[i], this->ExecTimesOnStep));

      if (referenceComparison) {
        this->learnReference(datasets[i]);
        refinementDetailsReference.push_back(
            make_pair(datasets[i], this->ExecTimesOnStepReference));
      }

      myFile << experimentHeaders[i] << csvSep << this->myTiming.timeComplete_;

      size_t lastStep = this->ExecTimesOnStep.size() - 1;
      myFile << csvSep << this->ExecTimesOnStep[lastStep].second;

      if (referenceComparison) {
        myFile << csvSep << this->referenceTiming.timeComplete_;
        myFile << csvSep << this->ExecTimesOnStepReference[lastStep].second;
        myFile << csvSep << this->referenceTiming.timeComplete_ / this->myTiming.timeComplete_;
        myFile << csvSep
               << this->ExecTimesOnStepReference[lastStep].second /
                      this->ExecTimesOnStep[lastStep].second;
      }

      myFile << std::endl;
    }

    myFile.close();

    this->writeRefinementResults(experimentDetailsFile, kernelMetaInformation, refinementDetails,
                                 refinementDetailsReference, referenceComparison);
  }
}

void MetaLearner::regularGridSpeedup(
    sgpp::datadriven::OperationMultipleEvalConfiguration &operationConfiguration,
    std::vector<size_t> dimList, std::vector<size_t> levelList, size_t instances,
    std::string metaInformation, std::string experimentName) {
  std::string kernelMetaInformation =
      metaInformation + " kernel: " + operationConfiguration.getName();

  std::ofstream myFile;
  std::string experimentFile =
      "results/data/" + experimentName + "_" + operationConfiguration.getName() + ".csv";
  remove(experimentFile.c_str());

  myFile.open(experimentFile, std::ios::out | std::ios::app);
  myFile << "# " << kernelMetaInformation << std::endl;

  myFile << "dim" << csvSep;
  myFile << "level" << csvSep;
  myFile << "duration" << csvSep;
  myFile << "durationRef" << csvSep;
  myFile << "speedup" << std::endl;

  for (size_t dim : dimList) {
    for (size_t level : levelList) {
      double duration;
      double durationReference;
      this->testRegular(operationConfiguration, dim, level, instances, duration, durationReference);

      myFile << dim << csvSep << level << csvSep << duration << csvSep << durationReference
             << csvSep;
      myFile << (durationReference / duration) << std::endl;
    }

    myFile << std::endl;
  }

  myFile.close();
}

void MetaLearner::appendToPerformanceRun(
    std::string fileName, std::string changingRowName, std::string currentValues,
    std::vector<sgpp::datadriven::OperationMultipleEvalConfiguration *> operationConfigurations,
    std::vector<std::string> datasets, std::vector<std::string> datasetNames,
    std::string metaInformation, bool removeOld) {
  if (removeOld) {
    std::ifstream fileForTest(fileName.c_str());

    if (!fileForTest.good()) {
      fileForTest.close();
      std::cout << "info: no old file to delete" << std::endl;
    } else {
      fileForTest.close();
      int returnCode = remove(fileName.c_str());

      if (returnCode == 0) {
        std::cout << "old result file deleted" << std::endl;
      } else {
        std::cout << "error: should remove old file, but couldn't, return code: " << returnCode
                  << std::endl;
        throw;
      }
    }

    std::ofstream myFile;
    myFile.open(fileName, std::ios::out | std::ios::app);
    myFile << "# " << metaInformation << std::endl;
    myFile << changingRowName;

    for (sgpp::datadriven::OperationMultipleEvalConfiguration *operationConfiguration :
         operationConfigurations) {
      for (std::string datasetName : datasetNames) {
        myFile << csvSep << operationConfiguration->getName() << "/" << datasetName;
      }
    }

    myFile << std::endl;
    myFile.close();
  }

  std::ofstream myFile;
  myFile.open(fileName, std::ios::out | std::ios::app);
  std::cout << "filename: " << fileName << std::endl;
  myFile << currentValues;

  for (sgpp::datadriven::OperationMultipleEvalConfiguration *operationConfiguration :
       operationConfigurations) {
    for (std::string dataset : datasets) {
      this->learn(*operationConfiguration, dataset);
      myFile << csvSep << this->myTiming.timeComplete_;
    }
  }

  myFile << std::endl;
  myFile.close();
}

double fRand(double fMin, double fMax) {
  double f = static_cast<double>(rand()) / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

void MetaLearner::testRegular(
    sgpp::datadriven::OperationMultipleEvalConfiguration &operationConfiguration, size_t dim,
    size_t level, size_t instances, double &duration, double &durationReference) {
  srand(static_cast<unsigned int>(time(NULL)));
  this->gridConfig.dim_ = dim;

  base::DataMatrix testTrainingData(instances, dim);

  for (size_t i = 0; i < instances; i++) {
    for (size_t d = 0; d < dim; d++) {
      // testTrainingData.set(i, d, unif(re));
      testTrainingData.set(i, d, fRand(0.0, 1.0));
    }
  }

  bool isRegression = true;
  LearnerLeastSquaresIdentity *learner =
      new LearnerLeastSquaresIdentity(isRegression, this->verbose);
  learner->setImplementation(operationConfiguration);

  LearnerLeastSquaresIdentity *learnerReference =
      new LearnerLeastSquaresIdentity(isRegression, this->verbose);
  sgpp::datadriven::OperationMultipleEvalConfiguration referenceOperationConfiguration(
      OperationMultipleEvalType::STREAMING, OperationMultipleEvalSubType::DEFAULT, "STREAMING");
  learnerReference->setImplementation(referenceOperationConfiguration);

  duration = learner->testRegular(this->gridConfig, testTrainingData);
  durationReference = learnerReference->testRegular(this->gridConfig, testTrainingData);
}

LearnerTiming MetaLearner::getLearnerTiming() { return this->myTiming; }

LearnerTiming MetaLearner::getLearnerReferenceTiming() { return this->referenceTiming; }
}  // namespace datadriven
}  // namespace sgpp
