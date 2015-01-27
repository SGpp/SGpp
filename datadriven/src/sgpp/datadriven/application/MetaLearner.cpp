// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include "MetaLearner.hpp"

#include <random>

#include <fstream>
#include <iomanip>
#include <cstdio>

#include "LearnerLeastSquaresIdentity.hpp"
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

using namespace SGPP::base;

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace datadriven {

    MetaLearner::MetaLearner(SGPP::solver::SLESolverConfiguration solverConfig,
                             SGPP::solver::SLESolverConfiguration solverFinalStep, SGPP::base::AdpativityConfiguration adaptivityConfiguration,
                             size_t baseLevel, double lambda, bool verbose) {
      this->csvSep = "& ";
      this->solverConfig = solverConfig;
      this->solverFinalStep = solverFinalStep;
      this->adaptivityConfiguration = adaptivityConfiguration;
      this->baseLevel = baseLevel;
      this->lambda = lambda;
      this->verbose = verbose;
      this->instances = 0;
      this->dim = 0;
    }

    void MetaLearner::learn(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
                            std::string datasetFileName) {

      Dataset dataset = ARFFTools::readARFF(datasetFileName);
      this->dim = dataset.getDimension();
      this->instances = dataset.getNumberInstances();

      if (verbose) {
        std::cout << "instances: " << this->instances << std::endl;
      }

      // Set grid config
      this->gridConfig.dim_ = this->dim;
      this->gridConfig.type_ = SGPP::base::Linear;
      this->gridConfig.level_ = static_cast<int>(this->baseLevel);

      DataVector* classesVector = dataset.getClasses();
      DataMatrix* trainingData = dataset.getTrainingData();

      bool isRegression = true; // treat everything as if it were a regression, as classification is not fully supported by Learner
      LearnerLeastSquaresIdentity* myLearner = new LearnerLeastSquaresIdentity(isRegression, this->verbose);
      myLearner->setImplementation(operationConfiguration);

      LearnerTiming timings = myLearner->train(*trainingData, *classesVector, this->gridConfig, this->solverConfig,
                              this->solverFinalStep, this->adaptivityConfiguration, false, this->lambda);

      this->myTiming = timings;
      this->ExecTimesOnStep = myLearner->getRefinementExecTimes();

      if (this->myLearner != nullptr) {
        delete this->myLearner;
      }

      this->myLearner = myLearner;
    }

    void MetaLearner::learnReference(std::string fileName) {

      Dataset dataset = ARFFTools::readARFF(fileName);
      this->dim = dataset.getDimension();
      this->instances = dataset.getNumberInstances();

      if (verbose) {
        std::cout << "instances: " << this->instances << std::endl;
      }

      // Set grid config
      this->gridConfig.dim_ = this->dim;
      this->gridConfig.type_ = SGPP::base::Linear;
      this->gridConfig.level_ = static_cast<int>(this->baseLevel);

      DataVector* classesVector = dataset.getClasses();
      DataMatrix* trainingData = dataset.getTrainingData();

      // treat everything as if it were a regression, as classification is not fully supported by Learner
      bool isRegression = true;
      LearnerLeastSquaresIdentity* referenceLearner = new LearnerLeastSquaresIdentity(isRegression, this->verbose);
      SGPP::datadriven::OperationMultipleEvalConfiguration operationConfiguration;
      operationConfiguration.type = OperationMultipleEvalType::DEFAULT;
      operationConfiguration.subType = OperationMultipleEvalSubType::DEFAULT;
      operationConfiguration.name = "STREAMING";
      referenceLearner->setImplementation(operationConfiguration);

      LearnerTiming timings = referenceLearner->train(*trainingData, *classesVector, gridConfig, solverConfig,
                              solverFinalStep, adaptivityConfiguration, false, lambda);
      this->referenceTiming = timings;
      this->ExecTimesOnStepReference = referenceLearner->getRefinementExecTimes();

      //referenceLearner->dumpFunction("referenceGridFunction", 50);
      if (this->referenceLearner != nullptr) {
        delete this->referenceLearner;
      }

      this->referenceLearner = referenceLearner;
    }

    //learn and test against test dataset and measure hits/mse
    void MetaLearner::learnAndTest(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
                                   std::string datasetFileName, std::string testFileName, bool isBinaryClassification) {

      //always to this first
      this->learn(operationConfiguration, datasetFileName);

      Dataset testDataset = ARFFTools::readARFF(testFileName);
      size_t testDim = testDataset.getDimension();
      size_t testInstances = testDataset.getNumberInstances();

      DataVector* testClassesVector = testDataset.getClasses();
      DataMatrix* testTrainingData = testDataset.getTrainingData();

      if (verbose && testDim != dim) {
        std::cout << "dim of test dataset and training dataset doesn't match" << std::endl;
      }

      if (verbose) {
        std::cout << "computing classes of test dataset" << std::endl;
      }

      DataVector computedClasses = this->myLearner->predict(*testTrainingData);

      if (verbose) {
        std::cout << "classes computed" << std::endl;
      }

      if (!isBinaryClassification) {
        double mse = 0.0;

        for (size_t i = 0; i < computedClasses.getSize(); i++) {
          double diff = testClassesVector->get(i) - computedClasses.get(i);
          mse += diff * diff;
        }

        mse = mse / (double) testInstances;

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

          if (classToken == testClassesVector->get(i)) {
            hits += 1;
          }
        }

        if (verbose) {
          std::cout << "hits (%): " << ((double) hits / (double) testInstances) << std::endl;
        }
      }
    }

    //learn and test against the streaming implemenation
    double MetaLearner::learnAndCompare(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
                                        std::string datasetFileName, size_t gridGranularity, double tolerance) {
      //always do this first
      this->learn(operationConfiguration, datasetFileName);
      this->learnReference(datasetFileName);

      DataMatrix testTrainingData(0, this->dim);

      double increment = 1.0 / static_cast<double>(gridGranularity);

      size_t index = 0;
      DataVector testPoint(dim);

      for (size_t i = 0; i < dim; i++) {
        testPoint[i] = 0.0 + increment;
      }

      testTrainingData.appendRow(testPoint);

      unsigned int testInstanceCounter = 0;

      while (index < dim) {
        // make 1.0 impossible
        double prior = testPoint[index] + increment;

        if (prior < 1.0) {

          testPoint[index] += increment;

          for (size_t i = 0; i < index; i++) {
            testPoint[i] = 0.0 + increment; //skip border
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
        std::cout << "testInstanceCounter: " << testInstanceCounter << std::endl;
        std::cout << "predicting..." << std::endl;
      }

      DataVector computedClasses = this->myLearner->predict(testTrainingData);

      if (verbose) {
        std::cout << "predicting... (reference)" << std::endl;
      }

      DataVector referenceClasses = this->referenceLearner->predict(testTrainingData);

      double squareSum = 0.0;

      for (size_t i = 0; i < computedClasses.getSize(); i++) {
        double temp = fabs(computedClasses.get(i) - referenceClasses.get(i));
        temp *= temp;
        squareSum += temp;

        if (verbose && fabs(computedClasses.get(i) - referenceClasses.get(i)) > tolerance) {
          std::cout << "computed: " << computedClasses.get(i) << " but reference is: " << referenceClasses.get(i)
                    << std::endl;
        }
      }

      squareSum = sqrt(squareSum);

      if (verbose) {
        std::cout << "sqrt: " << squareSum << std::endl;
      }

      return squareSum;
    }

    void MetaLearner::writeRefinementResults(std::string fileName, std::string fileHeader,
        std::vector<std::pair<std::string, std::vector<std::pair<size_t, double> > > > datasetDetails,
        std::vector<std::pair<std::string, std::vector<std::pair<size_t, double> > > > datasetDetailsReference,
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
      std::vector<SGPP::datadriven::OperationMultipleEvalConfiguration*> operationConfigurations,
      std::vector<std::string> datasets, std::vector<std::string> experimentHeaders, std::string metaInformation,
      std::string experimentName, bool referenceComparison) {
      for (SGPP::datadriven::OperationMultipleEvalConfiguration* operationConfiguration : operationConfigurations) {
        std::vector<std::pair<std::string, std::vector<std::pair<size_t, double> > > > refinementDetails;
        std::vector<std::pair<std::string, std::vector<std::pair<size_t, double> > > > refinementDetailsReference;

        std::string kernelMetaInformation = metaInformation + " kernel: " + operationConfiguration->name;

        std::ofstream myFile;
        std::string experimentFile = "results/data/overall_" + experimentName + "_" + operationConfiguration->name
                                     + ".csv";
        remove(experimentFile.c_str());
        std::string experimentDetailsFile = "results/data/refinement_" + experimentName + "_"
                                            + operationConfiguration->name + ".csv";
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
            refinementDetailsReference.push_back(make_pair(datasets[i], this->ExecTimesOnStepReference));
          }

          myFile << experimentHeaders[i] << csvSep << this->myTiming.timeComplete_;

          size_t lastStep = this->ExecTimesOnStep.size() - 1;
          myFile << csvSep << this->ExecTimesOnStep[lastStep].second;

          if (referenceComparison) {
            myFile << csvSep << this->referenceTiming.timeComplete_;
            myFile << csvSep << this->ExecTimesOnStepReference[lastStep].second;
            myFile << csvSep << this->referenceTiming.timeComplete_ / this->myTiming.timeComplete_;
            myFile << csvSep
                   << this->ExecTimesOnStepReference[lastStep].second / this->ExecTimesOnStep[lastStep].second;
          }

          myFile << std::endl;
        }

        myFile.close();

        this->writeRefinementResults(experimentDetailsFile, kernelMetaInformation, refinementDetails,
                                     refinementDetailsReference, referenceComparison);
      }
    }

    void MetaLearner::regularGridSpeedup(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
                                         std::vector<size_t> dimList, std::vector<size_t> levelList, size_t instances, std::string metaInformation,
                                         std::string experimentName) {

      std::string kernelMetaInformation = metaInformation + " kernel: " + operationConfiguration.name;

      std::ofstream myFile;
      std::string experimentFile = "results/data/" + experimentName + "_" + operationConfiguration.name + ".csv";
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

          myFile << dim << csvSep << level << csvSep << duration << csvSep << durationReference << csvSep;
          myFile << (durationReference / duration) << std::endl;
        }

        myFile << std::endl;
      }

      myFile.close();

    }

    void MetaLearner::appendToPerformanceRun(std::string fileName, std::string changingRowName, std::string currentValues,
        std::vector<SGPP::datadriven::OperationMultipleEvalConfiguration*> operationConfigurations,
        std::vector<std::string> datasets, std::vector<std::string> datasetNames, std::string metaInformation,
        bool removeOld) {

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
            std::cout << "error: should remove old file, but couldn't, return code: " << returnCode << std::endl;
            throw;
          }
        }

        std::ofstream myFile;
        myFile.open(fileName, std::ios::out | std::ios::app);
        myFile << "# " << metaInformation << std::endl;
        myFile << changingRowName;

        for (SGPP::datadriven::OperationMultipleEvalConfiguration* operationConfiguration : operationConfigurations) {
          for (std::string datasetName : datasetNames) {
            myFile << csvSep << operationConfiguration->name << "/" << datasetName;
          }
        }

        myFile << std::endl;
        myFile.close();
      }

      std::ofstream myFile;
      myFile.open(fileName, std::ios::out | std::ios::app);
      std::cout << "filename: " << fileName << std::endl;
      myFile << currentValues;

      for (SGPP::datadriven::OperationMultipleEvalConfiguration* operationConfiguration : operationConfigurations) {
        for (std::string dataset : datasets) {
          this->learn(*operationConfiguration, dataset);
          myFile << csvSep << this->myTiming.timeComplete_;
        }
      }

      myFile << std::endl;
      myFile.close();
    }

    double fRand(double fMin, double fMax) {
      double f = (double) rand() / RAND_MAX;
      return fMin + f * (fMax - fMin);
    }

    void MetaLearner::testRegular(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration, size_t dim,
                                  size_t level, size_t instances, double& duration, double& durationReference) {

      srand(static_cast<unsigned int>(time(NULL)));
      this->dim = dim;

      // Set grid config
      this->gridConfig.dim_ = dim;
      this->gridConfig.type_ = SGPP::base::Linear;
      this->gridConfig.level_ = static_cast<int>(level);

      DataMatrix testTrainingData(instances, dim);

      for (size_t i = 0; i < instances; i++) {
        for (size_t d = 0; d < dim; d++) {
          //testTrainingData.set(i, d, unif(re));
          testTrainingData.set(i, d, fRand(0.0, 1.0));
        }
      }

      bool isRegression = true;
      LearnerLeastSquaresIdentity* learner = new LearnerLeastSquaresIdentity(isRegression, this->verbose);
      learner->setImplementation(operationConfiguration);

      LearnerLeastSquaresIdentity* learnerReference = new LearnerLeastSquaresIdentity(isRegression,
          this->verbose);
      SGPP::datadriven::OperationMultipleEvalConfiguration referenceOperationConfiguration;
      referenceOperationConfiguration.type = OperationMultipleEvalType::STREAMING;
      referenceOperationConfiguration.subType = OperationMultipleEvalSubType::DEFAULT;
      referenceOperationConfiguration.name = "STREAMING";
      learnerReference->setImplementation(referenceOperationConfiguration);

      duration = learner->testRegular(this->gridConfig, testTrainingData);
      durationReference = learnerReference->testRegular(this->gridConfig, testTrainingData);
    }

  }
}

