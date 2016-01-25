// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <iostream>
#include <string>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/tools/TypesDatadriven.hpp>
#include <sgpp/datadriven/application/LearnerLeastSquaresIdentity.hpp>
#include <sgpp/solver/SLESolver.hpp>

#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace datadriven {

    class MetaLearner {
      private:

        size_t instances;

        std::string csvSep;

        bool verbose;

        std::shared_ptr<LearnerBase> myLearner;
        std::shared_ptr<LearnerBase> referenceLearner;

        base::RegularGridConfiguration gridConfig;

        solver::SLESolverConfiguration solverConfig;

        solver::SLESolverConfiguration solverFinalStep;

        base::AdpativityConfiguration adaptivityConfiguration;

        LearnerTiming myTiming;
        LearnerTiming referenceTiming;

        std::vector<std::pair<size_t, float_t> > ExecTimesOnStep;
        std::vector<std::pair<size_t, float_t> > ExecTimesOnStepReference;

        void writeRefinementResults(std::string fileName, std::string fileHeader,
                                    std::vector<std::pair<std::string, std::vector<std::pair<size_t, float_t> > > > datasetDetails,
                                    std::vector<std::pair<std::string, std::vector<std::pair<size_t, float_t> > > > datasetDetailsReference,
                                    bool referenceComparison);

        void optimizeLambdaLog_(SGPP::base::DataMatrix& dataset, SGPP::base::DataVector& datasetValues, size_t kFold,
                                size_t maxLevel, std::vector<base::DataMatrix>& trainingSets,
                                std::vector<base::DataVector>& trainingSetsValues, std::vector<base::DataMatrix>& testSets,
                                std::vector<base::DataVector>& testSetsValues, size_t curLevel, float_t lambdaLogStepSize,
                                float_t& bestLogLambda, float_t& bestMSE);

      public:
        MetaLearner() = delete;

        // gridConfig.dim is inferred from the dataset
        MetaLearner(SGPP::base::RegularGridConfiguration gridConfig, SGPP::solver::SLESolverConfiguration solverConfig,
                    SGPP::solver::SLESolverConfiguration solverFinalStep,
                    SGPP::base::AdpativityConfiguration adaptivityConfiguration, bool verbose = false);

        void learn(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
                   std::string& datasetFileName, float_t lambda, bool isRegression = true);

        void learnString(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration, std::string& content,
                         float_t lambda, bool isRegression = true);

        void learnReference(std::string& fileName, float_t lambda, bool isRegression = true);

        void learnReferenceString(std::string& content, float_t lambda, bool isRegression = true);

        //learn and test against test dataset and measure hits/mse
        void learnAndTest(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
                          std::string& datasetFileName, std::string& testFileName, bool isRegression = true);

        //learn and test against test dataset and measure hits/mse
        void learnAndTestString(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
                                std::string& dataContent, std::string& testContent, float_t lambda, bool isRegression = true);

        //learn and test against the streaming implementation
        float_t learnAndCompare(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
                                std::string& datasetFileName, float_t lambda, size_t gridGranularity);

        //learn and test against the streaming implementation
        float_t learnAndCompareString(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
                                      std::string& content, float_t lambda, size_t gridGranularity);

        std::shared_ptr<base::Grid> getLearnedGrid();

        LearnerTiming getLearnerTiming();

        LearnerTiming getLearnerReferenceTiming();

        float_t optimizeLambdaLog(SGPP::base::DataMatrix& dataset, SGPP::base::DataVector& values, size_t kFold,
                                  size_t maxLevel);

        void optimizeLambdaLog(SGPP::base::DataMatrix& dataset, SGPP::base::DataVector& values, size_t kFold,
                               size_t maxLevel, std::shared_ptr<base::Grid>& bestGrid, std::shared_ptr<base::DataVector>& bestAlpha,
                               float_t& lambdaOpt, datadriven::LearnerTiming& timing);

        void train(base::DataMatrix& train, base::DataVector& trainValues, float_t lambda,
                   std::shared_ptr<base::Grid>& grid, std::shared_ptr<base::DataVector>& alpha,
                   datadriven::LearnerTiming& timing);

        float_t calculateMSE(base::Grid& grid, base::DataVector& alpha, base::DataMatrix& testSubset,
                             base::DataVector& valuesTestSubset, bool verbose = false);

    };

  }
}
