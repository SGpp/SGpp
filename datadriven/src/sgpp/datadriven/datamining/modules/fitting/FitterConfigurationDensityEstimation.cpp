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
  FitterConfiguration::setupDefaults();
  // (Sebastian) Previously densityEstimationConfig was here set to
  // chol-decomp but has moved to the parent class FitterConfiguration
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
