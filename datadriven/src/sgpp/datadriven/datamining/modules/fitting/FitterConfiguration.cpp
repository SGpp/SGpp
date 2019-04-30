// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "FitterConfiguration.hpp"

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

const base::RegularGridConfiguration &FitterConfiguration::getGridConfig() const {
  return gridConfig;
}

const base::AdaptivityConfiguration &FitterConfiguration::getRefinementConfig() const {
  return adaptivityConfig;
}

const datadriven::CrossvalidationConfiguration &FitterConfiguration::getCrossvalidationConfig()
    const {
  return crossvalidationConfig;
}

const datadriven::DensityEstimationConfiguration &FitterConfiguration::getDensityEstimationConfig()
    const {
  return densityEstimationConfig;
}

const solver::SLESolverConfiguration &FitterConfiguration::getSolverRefineConfig() const {
  return solverRefineConfig;
}

const solver::SLESolverConfiguration &FitterConfiguration::getSolverFinalConfig() const {
  return solverFinalConfig;
}

const datadriven::RegularizationConfiguration &FitterConfiguration::getRegularizationConfig()
    const {
  return regularizationConfig;
}

const datadriven::OperationMultipleEvalConfiguration &FitterConfiguration::getMultipleEvalConfig()
    const {
  return multipleEvalConfig;
}

const datadriven::DatabaseConfiguration &FitterConfiguration::getDatabaseConfig() const {
  return databaseConfig;
}

const datadriven::LearnerConfiguration &FitterConfiguration::getLearnerConfig() const {
  return learnerConfig;
}

void FitterConfiguration::create(sgpp::datadriven::FitterConfiguration &m) {
    // Set Grid Config

    gridConfig.level_ = m.gridConfig.level_;
    gridConfig.boundaryLevel_ = m.gridConfig.boundaryLevel_;
    gridConfig.filename_ = m.gridConfig.filename_;
    gridConfig.generalType_ = m.gridConfig.generalType_;
    gridConfig.maxDegree_ = m.gridConfig.maxDegree_;
    gridConfig.dim_ = m.gridConfig.dim_;
    gridConfig.t_ = m.gridConfig.t_;
    gridConfig.type_ = m.gridConfig.type_;

    // Set Adaptivity Config

    adaptivityConfig.errorBasedRefinement = m.adaptivityConfig.errorBasedRefinement;
    adaptivityConfig.errorBufferSize = m.adaptivityConfig.errorBufferSize;
    adaptivityConfig.errorConvergenceThreshold = m.adaptivityConfig.errorConvergenceThreshold;
    adaptivityConfig.errorMinInterval = m.adaptivityConfig.errorMinInterval;
    adaptivityConfig.levelPenalize = m.adaptivityConfig.levelPenalize;
    adaptivityConfig.maxLevelType_ = m.adaptivityConfig.maxLevelType_;
    adaptivityConfig.noPoints_ = m.adaptivityConfig.noPoints_;
    adaptivityConfig.numRefinements_ = m.adaptivityConfig.numRefinements_;
    adaptivityConfig.percent_ = m.adaptivityConfig.percent_;
    adaptivityConfig.precomputeEvaluations = m.adaptivityConfig.precomputeEvaluations;
    adaptivityConfig.refinementFunctorType = m.adaptivityConfig.refinementFunctorType;
    adaptivityConfig.refinementPeriod = m.adaptivityConfig.refinementPeriod;
    adaptivityConfig.scalingCoefficients = m.adaptivityConfig.scalingCoefficients;
    adaptivityConfig.threshold_ = m.adaptivityConfig.threshold_;

    // Set Cross-Validation Config

    crossvalidationConfig.enable_ = m.crossvalidationConfig.enable_;
    crossvalidationConfig.kfold_ = m.crossvalidationConfig.kfold_;
    crossvalidationConfig.lambda_ = m.crossvalidationConfig.lambda_;
    crossvalidationConfig.lambdaEnd_ = m.crossvalidationConfig.lambdaEnd_;
    crossvalidationConfig.lambdaStart_ = m.crossvalidationConfig.lambdaStart_;
    crossvalidationConfig.logScale_ = m.crossvalidationConfig.logScale_;
    crossvalidationConfig.lambdaSteps_ = m.crossvalidationConfig.lambdaSteps_;
    crossvalidationConfig.seed_ = m.crossvalidationConfig.seed_;
    crossvalidationConfig.shuffle_ = m.crossvalidationConfig.shuffle_;
    crossvalidationConfig.silent_ = m.crossvalidationConfig.silent_;


    // Set Density Estimation Configuration

    densityEstimationConfig.type_ = m.densityEstimationConfig.type_;
    densityEstimationConfig.decomposition_ = m.densityEstimationConfig.decomposition_;
    densityEstimationConfig.iCholSweepsDecompose_ = m.densityEstimationConfig.iCholSweepsDecompose_;
    densityEstimationConfig.iCholSweepsRefine_ = m.densityEstimationConfig.iCholSweepsRefine_;
    densityEstimationConfig.iCholSweepsSolver_ = m.densityEstimationConfig.iCholSweepsSolver_;
    densityEstimationConfig.iCholSweepsUpdateLambda_ =
            m.densityEstimationConfig.iCholSweepsUpdateLambda_;


    // Set database Config

    databaseConfig.filepath = m.databaseConfig.filepath;

    // Set solver Refine Config


    solverRefineConfig.type_ = m.solverRefineConfig.type_;
    solverRefineConfig.threshold_ = m.solverRefineConfig.threshold_;
    solverRefineConfig.eps_ = m.solverRefineConfig.eps_;
    solverRefineConfig.maxIterations_ = m.solverRefineConfig.maxIterations_;
    solverRefineConfig.verbose_ = m.solverRefineConfig.verbose_;


    // Set Solver Final Config

    solverFinalConfig.verbose_ = m.solverFinalConfig.verbose_;
    solverFinalConfig.maxIterations_ = m.solverFinalConfig.maxIterations_;
    solverFinalConfig.eps_ = m.solverFinalConfig.eps_;
    solverFinalConfig.threshold_ = m.solverFinalConfig.threshold_;
    solverFinalConfig.type_ = m.solverFinalConfig.type_;


    // Set Regularization Config


    regularizationConfig.type_ = m.regularizationConfig.type_;
    regularizationConfig.lambda_ = m.regularizationConfig.lambda_;
    regularizationConfig.exponentBase_ = m.regularizationConfig.exponentBase_;
    regularizationConfig.l1Ratio_ = m.regularizationConfig.l1Ratio_;


    // Set multiple Eval  Config
    multipleEvalConfig = m.getMultipleEvalConfig();


    // Set learner Config

    learnerConfig.beta = m.learnerConfig.beta;
    learnerConfig.usePrior = m.learnerConfig.usePrior;
}

base::RegularGridConfiguration& FitterConfiguration::getGridConfig() {
  return const_cast<base::RegularGridConfiguration&>(
      static_cast<const FitterConfiguration&>(*this).getGridConfig());

const datadriven::ParallelConfiguration &FitterConfiguration::getParallelConfig() const {
  return parallelConfig;
}

base::RegularGridConfiguration &FitterConfiguration::getGridConfig() {
  return const_cast<base::RegularGridConfiguration &>(
      static_cast<const FitterConfiguration &>(*this).getGridConfig());
}

base::AdaptivityConfiguration &FitterConfiguration::getRefinementConfig() {
  return const_cast<base::AdaptivityConfiguration &>(
      static_cast<const FitterConfiguration &>(*this).getRefinementConfig());
}

datadriven::CrossvalidationConfiguration &FitterConfiguration::getCrossvalidationConfig() {
  return const_cast<datadriven::CrossvalidationConfiguration &>(
      static_cast<const FitterConfiguration &>(*this).getCrossvalidationConfig());
}

datadriven::DensityEstimationConfiguration &FitterConfiguration::getDensityEstimationConfig() {
  return const_cast<datadriven::DensityEstimationConfiguration &>(
      static_cast<const FitterConfiguration &>(*this).getDensityEstimationConfig());
}

solver::SLESolverConfiguration &FitterConfiguration::getSolverRefineConfig() {
  return const_cast<solver::SLESolverConfiguration &>(
      static_cast<const FitterConfiguration &>(*this).getSolverRefineConfig());
}

solver::SLESolverConfiguration &FitterConfiguration::getSolverFinalConfig() {
  return const_cast<solver::SLESolverConfiguration &>(
      static_cast<const FitterConfiguration &>(*this).getSolverFinalConfig());
}

datadriven::RegularizationConfiguration &FitterConfiguration::getRegularizationConfig() {
  return const_cast<datadriven::RegularizationConfiguration &>(
      static_cast<const FitterConfiguration &>(*this).getRegularizationConfig());
}

datadriven::OperationMultipleEvalConfiguration &FitterConfiguration::getMultipleEvalConfig() {
  return const_cast<datadriven::OperationMultipleEvalConfiguration &>(
      static_cast<const FitterConfiguration &>(*this).getMultipleEvalConfig());
}

void FitterConfiguration::setupDefaults() {
  gridConfig.type_ = sgpp::base::GridType::Linear;  // mirrors struct default
  gridConfig.dim_ = 0;
  gridConfig.level_ = 3;
  gridConfig.maxDegree_ = 1;  // mirrors struct default
  gridConfig.boundaryLevel_ = 0;  // mirrors struct default
  gridConfig.filename_ = "";
  gridConfig.t_ = 0.0;  // mirrors struct default

  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.threshold_ = 0.0;
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.noPoints_ = 0;
  adaptivityConfig.percent_ = 1.0;  // mirrors struct default
  adaptivityConfig.errorBasedRefinement = false;  // mirrors struct default
  adaptivityConfig.errorConvergenceThreshold = 0.001;  // mirrors struct default
  adaptivityConfig.errorBufferSize = 3;  // mirrors struct default
  adaptivityConfig.errorMinInterval = 0;  // mirrors struct default
  adaptivityConfig.refinementPeriod = 1;  // mirrors struct default
  adaptivityConfig.refinementFunctorType =
    sgpp::base::RefinementFunctorType::Surplus;  // mirrors struct default
  adaptivityConfig.precomputeEvaluations = true;  // mirrors struct default
  adaptivityConfig.levelPenalize = false;  // mirrors struct default
  adaptivityConfig.scalingCoefficients = std::vector<double>();  // mirrors struct default;

  crossvalidationConfig.enable_ = false;  // mirrors struct default
  crossvalidationConfig.kfold_ = 5;  // mirrors struct default
  crossvalidationConfig.seed_ = 0;
  crossvalidationConfig.shuffle_ = false;
  crossvalidationConfig.silent_ = false;
  crossvalidationConfig.lambda_ = 0.001;
  crossvalidationConfig.lambdaStart_ = 0.001;
  crossvalidationConfig.lambdaEnd_ = 0.001;
  crossvalidationConfig.lambdaSteps_ = 0;
  crossvalidationConfig.logScale_ = false;

  // (Sebastian) The following two values were previously set
  // in the subclass FitterConfigurationDensityEstimation but were moved here
  // to have all of the default value config in this file.
  densityEstimationConfig.type_ = sgpp::datadriven::DensityEstimationType::Decomposition;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::Chol;

  densityEstimationConfig.iCholSweepsDecompose_ = 4;  // mirrors struct default;
  densityEstimationConfig.iCholSweepsRefine_ = 4;  // mirrors struct default;
  densityEstimationConfig.iCholSweepsUpdateLambda_ = 2;  // mirrors struct default;
  densityEstimationConfig.iCholSweepsSolver_ = 2;  // mirrors struct default;

  databaseConfig.filepath = "";

  solverRefineConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverRefineConfig.eps_ = 1e-12;
  solverRefineConfig.maxIterations_ = 100;
  solverRefineConfig.threshold_ = 1e-12;
  solverRefineConfig.verbose_ = false;

  solverFinalConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverFinalConfig.eps_ = 1e-12;
  solverFinalConfig.maxIterations_ = 100;
  solverFinalConfig.threshold_ = 1e-12;
  solverFinalConfig.verbose_ = false;

  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 0.01;
  regularizationConfig.l1Ratio_ = 0.0;
  regularizationConfig.exponentBase_ = 1.0;

  learnerConfig.beta = 1.0;  // mirrors struct default
  learnerConfig.usePrior = false;  // mirrors struct default
}

}  // namespace datadriven
}  // namespace sgpp
