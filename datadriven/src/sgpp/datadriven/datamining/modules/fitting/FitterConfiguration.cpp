// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/tools/DMConfigTools.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

const base::GeneralGridConfiguration &FitterConfiguration::getGridConfig() const {
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
const datadriven::GeometryConfiguration &FitterConfiguration::getGeometryConfig() const {
  return geometryConfig;
}

const datadriven::ParallelConfiguration &FitterConfiguration::getParallelConfig() const {
  return parallelConfig;
}

base::GeneralGridConfiguration &FitterConfiguration::getGridConfig() {
  return const_cast<base::GeneralGridConfiguration &>(
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

void FitterConfiguration::dumpToStream(std::ostream &stream_out) const {
  stream_out << "\n~~~ gridConfig ~~~\n" << std::endl;
  DMConfigTools::dumpToStream(gridConfig);

  stream_out << "\n~~~ adaptivityConfig ~~~\n" << std::endl;
  DMConfigTools::dumpToStream(adaptivityConfig);

  stream_out << "\n~~~ crossvalidationConfig ~~~\n" << std::endl;
  DMConfigTools::dumpToStream(crossvalidationConfig);

  stream_out << "\n~~~ densityEstimationConfig ~~~\n" << std::endl;
  DMConfigTools::dumpToStream(densityEstimationConfig);

  stream_out << "\n~~~ databaseConfig ~~~\n" << std::endl;
  DMConfigTools::dumpToStream(databaseConfig);

  stream_out << "\n~~~ solverRefineConfig ~~~\n" << std::endl;
  DMConfigTools::dumpToStream(solverRefineConfig);

  stream_out << "\n~~~ solverFinalConfig ~~~\n" << std::endl;
  DMConfigTools::dumpToStream(solverFinalConfig);

  stream_out << "\n~~~ regularizationConfig ~~~\n" << std::endl;
  DMConfigTools::dumpToStream(regularizationConfig);

  stream_out << "\n~~~ learnerConfig ~~~\n" << std::endl;
  DMConfigTools::dumpToStream(learnerConfig);

  stream_out << "\n~~~ parallelConfig ~~~\n" << std::endl;
  DMConfigTools::dumpToStream(parallelConfig);

  // For now these are not included
  // stream_out << "\n~~~ geometryConfig ~~~\n" << std::endl;

  // stream_out << "\n~~~ multipleEvalConfig ~~~\n" << std::endl;
}

void FitterConfiguration::setupDefaults() {
  // (Sebastian Kreisel) The comments "mirrors struct default" are no longer applicable since all
  // structs now have default values that (should) match the ones set here. The comments are kept
  // for history / debugging.
  gridConfig.type_ = sgpp::base::GridType::Linear;  // mirrors struct default
  gridConfig.dim_ = 0;
  gridConfig.level_ = 3;
  gridConfig.levelVector_ = std::vector<size_t>();
  gridConfig.maxDegree_ = 1;      // mirrors struct default
  gridConfig.boundaryLevel_ = 0;  // mirrors struct default
  gridConfig.filename_ = "";
  gridConfig.t_ = 0.0;  // mirrors struct default

  adaptivityConfig.numRefinements_ = 1;
  adaptivityConfig.refinementThreshold_ = 0.0;
  adaptivityConfig.coarseningThreshold_ = 0.0;
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.numRefinementPoints_ = 1;
  adaptivityConfig.numCoarseningPoints_ = 1;
  adaptivityConfig.coarsenInitialPoints_ = false;
  adaptivityConfig.percent_ = 1.0;                      // mirrors struct default
  adaptivityConfig.errorBasedRefinement_ = false;       // mirrors struct default
  adaptivityConfig.errorConvergenceThreshold_ = 0.001;  // mirrors struct default
  adaptivityConfig.errorBufferSize_ = 3;                // mirrors struct default
  adaptivityConfig.errorMinInterval_ = 0;               // mirrors struct default
  adaptivityConfig.refinementPeriod_ = 1;               // mirrors struct default
  adaptivityConfig.refinementFunctorType_ =
      sgpp::base::RefinementFunctorType::Surplus;  // mirrors struct default
  adaptivityConfig.coarseningFunctorType_ =        // mirrors struct default
      sgpp::base::CoarseningFunctorType::Surplus;
  adaptivityConfig.precomputeEvaluations_ = true;                 // mirrors struct default
  adaptivityConfig.levelPenalize_ = false;                        // mirrors struct default
  adaptivityConfig.scalingCoefficients_ = std::vector<double>();  // mirrors struct default;

  crossvalidationConfig.enable_ = false;  // mirrors struct default
  crossvalidationConfig.kfold_ = 5;       // mirrors struct default
  crossvalidationConfig.seed_ = 0;
  crossvalidationConfig.shuffle_ = false;
  crossvalidationConfig.silent_ = false;
  crossvalidationConfig.lambda_ = 0.001;
  crossvalidationConfig.lambdaStart_ = 0.001;
  crossvalidationConfig.lambdaEnd_ = 0.001;
  crossvalidationConfig.lambdaSteps_ = 0;
  crossvalidationConfig.logScale_ = false;

  // (Sebastian) The following two values were previously set in the subclass
  // FitterConfigurationDensityEstimation but were moved here to have all of the default value
  // config in this file.
  densityEstimationConfig.type_ = sgpp::datadriven::DensityEstimationType::Decomposition;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;
  // Offline permutation is used per default
  densityEstimationConfig.useOfflinePermutation_ = true;

  densityEstimationConfig.iCholSweepsDecompose_ = 4;     // mirrors struct default
  densityEstimationConfig.iCholSweepsRefine_ = 4;        // mirrors struct default
  densityEstimationConfig.iCholSweepsUpdateLambda_ = 2;  // mirrors struct default
  densityEstimationConfig.iCholSweepsSolver_ = 2;        // mirrors struct default

  databaseConfig.filePath_ = "";

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
  regularizationConfig.lamda_start_ = 0.01;
  regularizationConfig.lambda_end_ = 0.01;
  regularizationConfig.lambda_steps_ = 0;
  regularizationConfig.lambda_log_scale_ = false;
  regularizationConfig.l1Ratio_ = 0.0;
  regularizationConfig.exponentBase_ = 1.0;
  regularizationConfig.optimizeLambda_ = false;
  regularizationConfig.optimizerTolerance_ = 1e-15;
  regularizationConfig.convergenceThreshold_ = 1e-5;
  regularizationConfig.intervalA_ = 1e-15;
  regularizationConfig.intervalB_ = 1.0;

  learnerConfig.learningRate_ = 1.0;  // mirrors struct default
  learnerConfig.usePrior_ = false;    // mirrors struct default

  // configure geometry configuration
  geometryConfig.dim_ = std::vector<std::vector<int64_t>>();
  geometryConfig.stencils_ = std::vector<sgpp::datadriven::StencilConfiguration>();
}

}  // namespace datadriven
}  // namespace sgpp
