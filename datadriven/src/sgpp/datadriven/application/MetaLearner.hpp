// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of
// distribution and use, please see the copyright notice
// provided with SG++ or at sgpp.sparsegrids.org.

#pragma once

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "sgpp/base/datatypes/DataVector.hpp"
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/grid/GridStorage.hpp"
#include "sgpp/base/grid/generation/GridGenerator.hpp"
#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/datadriven/tools/TypesDatadriven.hpp"
#include "sgpp/datadriven/application/LearnerLeastSquaresIdentity.hpp"
#include "sgpp/solver/SLESolver.hpp"
#include "sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp"
#include "sgpp/globaldef.hpp"

namespace SGPP {
namespace datadriven {

class MetaLearner {
 private:
  size_t instances;
  float_t lambda;

  std::string csvSep;

  bool verbose;

  LearnerBase* myLearner = nullptr;
  LearnerBase* referenceLearner = nullptr;

  SGPP::base::RegularGridConfiguration gridConfig;
  SGPP::solver::SLESolverConfiguration solverConfig;
  SGPP::solver::SLESolverConfiguration solverFinalStep;
  SGPP::base::AdpativityConfiguration adaptivityConfiguration;

  LearnerTiming myTiming;
  LearnerTiming referenceTiming;

  std::vector<std::pair<size_t, float_t> > ExecTimesOnStep;
  std::vector<std::pair<size_t, float_t> > ExecTimesOnStepReference;

  void writeRefinementResults(
      std::string fileName, std::string fileHeader,
      std::vector<std::pair<std::string, std::vector<std::pair<size_t, float_t> > > >
          datasetDetails,
      std::vector<std::pair<std::string, std::vector<std::pair<size_t, float_t> > > >
          datasetDetailsReference,
      bool referenceComparison);

 public:
  MetaLearner() = delete;

  // gridConfig.dim is inferred from the dataset
  MetaLearner(SGPP::base::RegularGridConfiguration gridConfig,
              SGPP::solver::SLESolverConfiguration solverConfig,
              SGPP::solver::SLESolverConfiguration solverFinalStep,
              SGPP::base::AdpativityConfiguration adaptivityConfiguration, float_t lambda,
              bool verbose = false);

  ~MetaLearner() {
    if (this->myLearner != nullptr) {
      delete this->myLearner;
    }

    if (this->referenceLearner != nullptr) {
      delete this->referenceLearner;
    }
  }

  void learn(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
             std::string& datasetFileName, bool isRegression = true);

  void learnString(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
                   std::string& content, bool isRegression = true);

  void learnReference(std::string& fileName, bool isRegression = true);

  void learnReferenceString(std::string& content, bool isRegression = true);

  // learn and test against test dataset and measure hits/mse
  void learnAndTest(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
                    std::string& datasetFileName, std::string& testFileName,
                    bool isRegression = true);

  // learn and test against test dataset and measure hits/mse
  void learnAndTestString(
      SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
      std::string& dataContent, std::string& testContent, bool isRegression = true);

  // learn and test against the streaming implementation
  float_t learnAndCompare(
      SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
      std::string& datasetFileName, size_t gridGranularity);

  // learn and test against the streaming implementation
  float_t learnAndCompareString(
      SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
      std::string& content, size_t gridGranularity);

  void refinementAndOverallPerformance(
      std::vector<SGPP::datadriven::OperationMultipleEvalConfiguration*> operationConfigurations,
      std::vector<std::string> datasets, std::vector<std::string> experimentHeaders,
      std::string metaInformation, std::string fileName, bool referenceComparison = false);

  void regularGridSpeedup(
      SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
      std::vector<size_t> dimList, std::vector<size_t> levelList, size_t instances,
      std::string metaInformation, std::string experimentName);

  void appendToPerformanceRun(
      std::string fileName, std::string changingRowName, std::string currentValues,
      std::vector<SGPP::datadriven::OperationMultipleEvalConfiguration*> operationConfigurations,
      std::vector<std::string> datasets, std::vector<std::string> datasetNames,
      std::string metaInformation, bool removeOld);

  void testRegular(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
                   size_t dim, size_t level, size_t instances, float_t& duration,
                   float_t& durationReference);

  SGPP::base::Grid& getLearnedGrid();

  LearnerTiming getLearnerTiming();

  LearnerTiming getLearnerReferenceTiming();
};
}  // namespace datadriven
}  // namespace SGPP
