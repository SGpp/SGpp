// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/AdaptivityThresholdTypeParser.hpp>
#include <sgpp/base/grid/CoarseningFunctorTypeParser.hpp>
#include <sgpp/base/grid/GeneralGridTypeParser.hpp>
#include <sgpp/base/grid/GridTypeParser.hpp>
#include <sgpp/base/grid/RefinementFunctorTypeParser.hpp>

#include <sgpp/datadriven/configuration/DensityEstimationTypeParser.hpp>
#include <sgpp/datadriven/configuration/MatrixDecompositionTypeParser.hpp>
#include <sgpp/datadriven/configuration/RegularizationTypeParser.hpp>

#include <sgpp/solver/SLESolverTypeParser.hpp>

namespace sgpp {
namespace datadriven {

class DMConfigTools {
 public:
  static void dumpToStream(const base::GeneralGridConfiguration& gridConfig,
                           std::ostream& stream_out = std::cout) {
    stream_out << "type \t\t\t" << base::GeneralGridTypeParser::toString(gridConfig.generalType_)
               << std::endl;
    stream_out << "type \t\t\t" << base::GridTypeParser::toString(gridConfig.type_) << std::endl;
    stream_out << "dim \t\t\t" << gridConfig.dim_ << std::endl;
    stream_out << "level \t\t\t" << gridConfig.level_ << std::endl;
    stream_out << "levelVector \t\t";
    for (auto i = gridConfig.levelVector_.begin(); i != gridConfig.levelVector_.end(); ++i)
      stream_out << *i << ' ';
    stream_out << std::endl;
    stream_out << "maxDegree \t\t" << gridConfig.maxDegree_ << std::endl;
    stream_out << "boundaryLevel \t\t" << gridConfig.boundaryLevel_ << std::endl;
    stream_out << "filename \t\t" << gridConfig.filename_ << std::endl;
    stream_out << "t \t\t\t" << gridConfig.t_ << std::endl;
  }

  static void dumpToStream(const base::RegularGridConfiguration& gridConfig,
                           std::ostream& stream_out = std::cout) {
    // General grid config
    stream_out << "type \t\t\t" << base::GeneralGridTypeParser::toString(gridConfig.generalType_)
               << std::endl;
    stream_out << "type \t\t\t" << base::GridTypeParser::toString(gridConfig.type_) << std::endl;
    stream_out << "dim \t\t\t" << gridConfig.dim_ << std::endl;
    stream_out << "level \t\t\t" << gridConfig.level_ << std::endl;
    stream_out << "levelVector \t\t";
    for (auto i = gridConfig.levelVector_.begin(); i != gridConfig.levelVector_.end(); ++i)
      stream_out << *i << ' ';
    stream_out << std::endl;
    stream_out << "maxDegree \t\t" << gridConfig.maxDegree_ << std::endl;
    stream_out << "boundaryLevel \t\t" << gridConfig.boundaryLevel_ << std::endl;
    stream_out << "filename \t\t" << gridConfig.filename_ << std::endl;
    stream_out << "t \t\t\t" << gridConfig.t_ << std::endl;
  }

  static void dumpToStream(const base::CombiGridConfiguration& gridConfig,
                           std::ostream& stream_out = std::cout) {
    // General grid config
    stream_out << "type \t\t\t" << base::GeneralGridTypeParser::toString(gridConfig.generalType_)
               << std::endl;
    stream_out << "type \t\t\t" << base::GridTypeParser::toString(gridConfig.type_) << std::endl;
    stream_out << "dim \t\t\t" << gridConfig.dim_ << std::endl;
    stream_out << "level \t\t\t" << gridConfig.level_ << std::endl;
    stream_out << "levelVector \t\t";
    for (auto i = gridConfig.levelVector_.begin(); i != gridConfig.levelVector_.end(); ++i)
      stream_out << *i << ' ';
    stream_out << std::endl;
    stream_out << "maxDegree \t\t" << gridConfig.maxDegree_ << std::endl;
    stream_out << "boundaryLevel \t\t" << gridConfig.boundaryLevel_ << std::endl;
    stream_out << "filename \t\t" << gridConfig.filename_ << std::endl;
    stream_out << "t \t\t\t" << gridConfig.t_ << std::endl;

    // Additional config
    stream_out << "levels \t\t\t";
    for (auto i = gridConfig.levels.begin(); i != gridConfig.levels.end(); ++i)
      stream_out << *i << ' ';
    stream_out << std::endl;
  }

  static void dumpToStream(const base::AdaptivityConfiguration& adaptivityConfig,
                           std::ostream& stream_out = std::cout) {
    stream_out << "numRefinements \t\t\t" << adaptivityConfig.numRefinements_ << std::endl;
    stream_out << "thresholdType \t\t\t"
               << base::AdaptivityThresholdTypeParser::toString(adaptivityConfig.thresholdType_)
               << std::endl;
    stream_out << "refinementThreshold \t\t" << adaptivityConfig.refinementThreshold_ << std::endl;
    stream_out << "refinementThreshold \t\t" << adaptivityConfig.coarseningThreshold_ << std::endl;
    stream_out << std::boolalpha << "coarsenInitialPoints \t\t"
               << adaptivityConfig.coarsenInitialPoints_ << std::endl;
    stream_out << std::boolalpha << "maxLevelType \t\t\t" << adaptivityConfig.maxLevelType_
               << std::endl;
    stream_out << "numRefinementPoints \t\t" << adaptivityConfig.numRefinementPoints_ << std::endl;
    stream_out << "numCoarseningPoints \t\t" << adaptivityConfig.numCoarseningPoints_ << std::endl;
    stream_out << "percent \t\t\t" << adaptivityConfig.percent_ << std::endl;
    stream_out << std::boolalpha << "errorBasedRefinement \t\t"
               << adaptivityConfig.errorBasedRefinement_ << std::endl;
    stream_out << "errorConvergenceThreshold \t" << adaptivityConfig.errorConvergenceThreshold_
               << std::endl;
    stream_out << "errorBufferSize \t\t" << adaptivityConfig.errorBufferSize_ << std::endl;
    stream_out << "errorMinInterval \t\t" << adaptivityConfig.errorMinInterval_ << std::endl;
    stream_out << "refinementPeriod \t\t" << adaptivityConfig.refinementPeriod_ << std::endl;
    stream_out << "refinementFunctorType \t\t" << base::RefinementFunctorTypeParser::toString(
                                                      adaptivityConfig.refinementFunctorType_)
               << std::endl;
    stream_out << "coarseningFunctorType \t\t" << base::CoarseningFunctorTypeParser::toString(
                                                      adaptivityConfig.coarseningFunctorType_)
               << std::endl;
    stream_out << std::boolalpha << "precomputeEvaluations \t\t"
               << adaptivityConfig.precomputeEvaluations_ << std::endl;
    stream_out << std::boolalpha << "levelPenalize \t\t\t" << adaptivityConfig.levelPenalize_
               << std::endl;
    stream_out << "scalingCoefficients \t\t";
    for (auto i = adaptivityConfig.scalingCoefficients_.begin();
         i != adaptivityConfig.scalingCoefficients_.end(); ++i)
      stream_out << *i << ' ';
    stream_out << std::endl;
  }

  static void dumpToStream(const CrossvalidationConfiguration& crossvalidationConfig,
                           std::ostream& stream_out = std::cout) {
    stream_out << "enable \t\t\t" << std::boolalpha << crossvalidationConfig.enable_ << std::endl;
    stream_out << "kfold \t\t\t" << crossvalidationConfig.kfold_ << std::endl;
    stream_out << "seed \t\t\t" << crossvalidationConfig.seed_ << std::endl;
    stream_out << "shuffle \t\t" << std::boolalpha << crossvalidationConfig.shuffle_ << std::endl;
    stream_out << "silent \t\t\t" << std::boolalpha << crossvalidationConfig.silent_ << std::endl;
    stream_out << "lambda \t\t\t" << crossvalidationConfig.lambda_ << std::endl;
    stream_out << "lambdaEnd \t\t" << crossvalidationConfig.lambdaEnd_ << std::endl;
    stream_out << "lambdaStart \t\t" << crossvalidationConfig.lambdaStart_ << std::endl;
    stream_out << "lambdaSteps \t\t" << crossvalidationConfig.lambdaSteps_ << std::endl;
    stream_out << "logScale \t\t" << std::boolalpha << crossvalidationConfig.logScale_ << std::endl;
  }

  static void dumpToStream(const DensityEstimationConfiguration& densityEstimationConfig,
                           std::ostream& stream_out = std::cout) {
    stream_out << "type \t\t\t\t"
               << DensityEstimationTypeParser::toString(densityEstimationConfig.type_) << std::endl;
    stream_out << "decomposition \t\t\t" << datadriven::MatrixDecompositionTypeParser::toString(
                                                densityEstimationConfig.decomposition_)
               << std::endl;
    stream_out << "useOfflinePermutation \t\t" << std::boolalpha
               << densityEstimationConfig.useOfflinePermutation_ << std::endl;
    stream_out << "normalize \t\t\t" << std::boolalpha << densityEstimationConfig.normalize_
               << std::endl;
    stream_out << "iCholSweepsDecompose \t\t" << densityEstimationConfig.iCholSweepsDecompose_
               << std::endl;
    stream_out << "iCholSweepsRefine \t\t" << densityEstimationConfig.iCholSweepsRefine_
               << std::endl;
    stream_out << "iCholSweepsUpdateLambda \t" << densityEstimationConfig.iCholSweepsUpdateLambda_
               << std::endl;
    stream_out << "iCholSweepsSolver \t\t" << densityEstimationConfig.iCholSweepsSolver_
               << std::endl;
  }

  static void dumpToStream(const DatabaseConfiguration& databaseConfig,
                           std::ostream& stream_out = std::cout) {
    stream_out << "filePath \t\t" << databaseConfig.filePath_ << std::endl;
  }

  static void dumpToStream(const solver::SLESolverConfiguration& solverConfig,
                           std::ostream& stream_out = std::cout) {
    stream_out << "type \t\t\t" << solver::SLESolverTypeParser::toString(solverConfig.type_)
               << std::endl;
    stream_out << "eps \t\t\t" << solverConfig.eps_ << std::endl;
    stream_out << "maxIterations \t\t" << solverConfig.maxIterations_ << std::endl;
    stream_out << "threshold \t\t" << solverConfig.threshold_ << std::endl;
    stream_out << "verbose \t\t" << std::boolalpha << solverConfig.verbose_ << std::endl;
  }

  static void dumpToStream(const RegularizationConfiguration& regularizationConfig,
                           std::ostream& stream_out = std::cout) {
    stream_out << "type \t\t\t" << RegularizationTypeParser::toString(regularizationConfig.type_)
               << std::endl;
    stream_out << "lambda \t\t\t" << regularizationConfig.lambda_ << std::endl;
    stream_out << "l1Ratio \t\t" << regularizationConfig.l1Ratio_ << std::endl;
    stream_out << "exponentBase \t\t" << regularizationConfig.exponentBase_ << std::endl;
    stream_out << "lambda_start \t\t" << regularizationConfig.lamda_start_ << std::endl;
    stream_out << "lambda_end \t\t" << regularizationConfig.lambda_end_ << std::endl;
    stream_out << "lambda_steps \t\t" << regularizationConfig.lambda_steps_ << std::endl;
    stream_out << "lambda_log_scale \t" << std::boolalpha << regularizationConfig.lambda_log_scale_
               << std::endl;
    stream_out << "optimizeLambda \t\t" << std::boolalpha << regularizationConfig.optimizeLambda_
               << std::endl;
    stream_out << "optimizerTolerance \t" << regularizationConfig.optimizerTolerance_ << std::endl;
    stream_out << "convergenceThreshold \t" << regularizationConfig.convergenceThreshold_
               << std::endl;
    stream_out << "intervalA \t\t" << regularizationConfig.intervalA_ << std::endl;
    stream_out << "intervalB \t\t" << regularizationConfig.intervalB_ << std::endl;
  }

  static void dumpToStream(const LearnerConfiguration& learnerConfig,
                           std::ostream& stream_out = std::cout) {
    stream_out << "learningRate \t\t" << learnerConfig.learningRate_ << std::endl;
    stream_out << "usePrior \t\t" << std::boolalpha << learnerConfig.usePrior_ << std::endl;
  }

  static void dumpToStream(const ParallelConfiguration& parallelConfig,
                           std::ostream& stream_out = std::cout) {
    stream_out << "scalapackEnabled \t" << std::boolalpha << parallelConfig.scalapackEnabled_
               << std::endl;
    stream_out << "processRows \t\t" << parallelConfig.processRows_ << std::endl;
    stream_out << "processCols \t\t" << parallelConfig.processCols_ << std::endl;
    stream_out << "rowBlockSize \t\t" << parallelConfig.rowBlockSize_ << std::endl;
    stream_out << "columnBlockSize \t" << parallelConfig.columnBlockSize_ << std::endl;
  }
};

} /* namespace datadriven */
} /* namespace sgpp */
