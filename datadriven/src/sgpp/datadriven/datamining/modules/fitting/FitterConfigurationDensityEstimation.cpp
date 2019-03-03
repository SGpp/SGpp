/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FitterConfigurationDensityEstimation.cpp
 *
 *  Created on: Jan 02, 2018
 *      Author: Kilian Röhner
 */

#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/solver/TypesSolver.hpp>

namespace sgpp {
namespace datadriven {

FitterConfiguration *FitterConfigurationDensityEstimation::clone() const {
  return new FitterConfigurationDensityEstimation(*this);
}

void FitterConfigurationDensityEstimation::setupDefaults() {
  // configure initial grid
  gridConfig.dim_ = 0;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  // configure adaptive refinement
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.refinementPeriod = 1;
  adaptivityConfig.errorBasedRefinement = false;

  // configure regularization
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 0.01;

  // configure crossvalidation
  crossvalidationConfig.enable_ = false;

  // configure density estimation
  densityEstimationConfig.type_ = sgpp::datadriven::DensityEstimationType::Decomposition;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::Chol;

  // intialize empty database
  databaseConfig.filepath = "";

  // configure learner
  learnerConfig.beta = 1.0;

  // Configure solvers for CG
  solverFinalConfig.maxIterations_ = 100;
  solverFinalConfig.eps_ = 1e-14;
  solverFinalConfig.threshold_ = 1e-14;
  solverFinalConfig.type_ = sgpp::solver::SLESolverType::CG;

  solverRefineConfig.maxIterations_ = 100;
  solverRefineConfig.eps_ = 1e-14;
  solverRefineConfig.threshold_ = 1e-14;
  solverRefineConfig.type_ = sgpp::solver::SLESolverType::CG;
}

void FitterConfigurationDensityEstimation::readParams(const DataMiningConfigParser &parser) {
  setupDefaults();

  parser.getFitterGridConfig(gridConfig, gridConfig);
  //parser.getFitterAdaptivityConfig(adaptivityConfig, adaptivityConfig);
  parser.getFitterCoarseningConfig(coarseningConfig, coarseningConfig);
  parser.getFitterSolverRefineConfig(solverRefineConfig, solverRefineConfig);
  parser.getFitterSolverFinalConfig(solverFinalConfig, solverFinalConfig);
  parser.getFitterRegularizationConfig(regularizationConfig, regularizationConfig);
  parser.getFitterDensityEstimationConfig(densityEstimationConfig, densityEstimationConfig);
  parser.getFitterDatabaseConfig(databaseConfig, databaseConfig);
  parser.getFitterLearnerConfig(learnerConfig, learnerConfig);
}
} /* namespace datadriven */
} /* namespace sgpp */
