// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "FitterConfiguration.hpp"

#include <string>

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

void FitterConfiguration::dumpParamsToStdout(std::ostream &stream_out) const {
  stream_out << "\n~~~ gridConfig ~~~\n" << std::endl;

  stream_out << "boundaryLevel \t\t" << gridConfig.boundaryLevel_ << std::endl;
  stream_out << "dim \t\t\t" << gridConfig.dim_ << std::endl;
  stream_out << "filename \t\t" << gridConfig.filename_ << std::endl;
  // stream_out << gridConfig.generalType_ << std::endl;
  stream_out << "level \t\t\t" << gridConfig.level_ << std::endl;
  stream_out << "maxDegree \t\t" << gridConfig.maxDegree_ << std::endl;
  stream_out << "t \t\t\t" << gridConfig.t_ << std::endl;
  // stream_out << gridConfig.type_ << std::endl;

  stream_out << "\n~~~ adaptivityConfig ~~~\n" << std::endl;

  stream_out << "maxLevelType \t\t" << adaptivityConfig.maxLevelType_ << std::endl;
  stream_out << "noPoints \t\t" << adaptivityConfig.noPoints_ << std::endl;
  stream_out << "numRefinements \t\t" << adaptivityConfig.numRefinements_ << std::endl;
  stream_out << "percent \t\t" << adaptivityConfig.percent_ << std::endl;

  stream_out << "\n~~~ crossvalidationConfig ~~~\n" << std::endl;

  stream_out << "enable \t\t\t" << crossvalidationConfig.enable_ << std::endl;
  stream_out << "kfold \t\t\t" << crossvalidationConfig.kfold_ << std::endl;
  stream_out << "lambdaEnd \t\t" << crossvalidationConfig.lambdaEnd_ << std::endl;
  stream_out << "lambdaStart \t\t" << crossvalidationConfig.lambdaStart_ << std::endl;
  stream_out << "lambdaSteps \t\t" << crossvalidationConfig.lambdaSteps_ << std::endl;
  stream_out << "lambda \t\t\t" << crossvalidationConfig.lambda_ << std::endl;
  stream_out << "logScale \t\t" << crossvalidationConfig.logScale_ << std::endl;
  stream_out << "seed \t\t\t" << crossvalidationConfig.seed_ << std::endl;
  stream_out << "shuffle \t\t" << crossvalidationConfig.shuffle_ << std::endl;
  stream_out << "silent \t\t\t" << crossvalidationConfig.silent_ << std::endl;

  stream_out << "\n~~~ densityEstimationConfig ~~~\n" << std::endl;

  // stream_out << densityEstimationConfig.decomposition_ << std::endl;
  stream_out << "iCholSweepsDecompose \t" << densityEstimationConfig.iCholSweepsDecompose_
             << std::endl;
  stream_out << "iCholSweepsRefine \t" << densityEstimationConfig.iCholSweepsRefine_ << std::endl;
  stream_out << "iCholSweepsSolver \t" << densityEstimationConfig.iCholSweepsSolver_ << std::endl;
  stream_out << "iCholSweepsUpdateLambda " << densityEstimationConfig.iCholSweepsUpdateLambda_
             << std::endl;
  // stream_out << densityEstimationConfig.type_ << std::endl;

  stream_out << "\n~~~ databaseConfig ~~~\n" << std::endl;

  stream_out << "filepath \t\t" << databaseConfig.filepath << std::endl;

  stream_out << "\n~~~ solverFinalConfig ~~~\n" << std::endl;

  stream_out << "eps \t\t\t" << solverFinalConfig.eps_ << std::endl;
  stream_out << "maxIterations \t" << solverFinalConfig.maxIterations_ << std::endl;
  stream_out << "threshold \t\t" << solverFinalConfig.threshold_ << std::endl;
  // stream_out << solverFinalConfig.type_ << std::endl;
  stream_out << "verbose \t\t" << solverFinalConfig.verbose_ << std::endl;

  stream_out << "\n~~~ solverRefineConfig ~~~\n" << std::endl;

  stream_out << "eps \t\t\t" << solverRefineConfig.eps_ << std::endl;
  stream_out << "maxIterations \t\t" << solverRefineConfig.maxIterations_ << std::endl;
  stream_out << "threshold \t\t" << solverRefineConfig.threshold_ << std::endl;
  // stream_out << solverRefineConfig.type_ << std::endl;
  stream_out << "verbose \t\t" << solverRefineConfig.verbose_ << std::endl;

  stream_out << "\n~~~ regularizationConfig ~~~\n" << std::endl;

  stream_out << "exponentBase \t\t" << regularizationConfig.exponentBase_ << std::endl;
  stream_out << "l1Ratio \t\t" << regularizationConfig.l1Ratio_ << std::endl;
  stream_out << "lambda \t\t\t" << regularizationConfig.lambda_ << std::endl;
  // stream_out << regularizationConfig.type_ << std::endl;

  stream_out << "\n~~~ multipleEvalConfig ~~~\n" << std::endl;

  // stream_out << multipleEvalConfig.getMPIType() << std::endl;
  stream_out << "getName() \t\t" << multipleEvalConfig.getName() << std::endl;

  stream_out << "\n~~~ learnerConfig ~~~\n" << std::endl;

  stream_out << "beta \t\t\t" << learnerConfig.beta << std::endl;
  stream_out << "usePrior \t\t" << learnerConfig.usePrior << std::endl;
}

}  // namespace datadriven
}  // namespace sgpp
