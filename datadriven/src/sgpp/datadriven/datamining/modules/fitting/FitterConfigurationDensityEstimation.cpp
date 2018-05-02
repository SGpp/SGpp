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

namespace sgpp {
namespace datadriven {

FitterConfiguration* FitterConfigurationDensityEstimation::clone() const {
  return new FitterConfigurationDensityEstimation(*this);
}

void FitterConfigurationDensityEstimation::setupDefaults() {
  // configure initial grid
  gridConfig.dim_ = 0;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  // configure adaptive refinement
  adaptivityConfig.numRefinements_ = 0;

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
}

void FitterConfigurationDensityEstimation::readParams(const DataMiningConfigParser& parser) {
  setupDefaults();

  parser.getFitterGridConfig(gridConfig, gridConfig);
  parser.getFitterAdaptivityConfig(adaptivityConfig, adaptivityConfig);
  parser.getFitterSolverRefineConfig(solverRefineConfig, solverRefineConfig);
  parser.getFitterSolverFinalConfig(solverFinalConfig, solverFinalConfig);
  parser.getFitterRegularizationConfig(regularizationConfig, regularizationConfig);
  parser.getFitterDatabaseConfig(databaseConfig, databaseConfig);
}
} /* namespace datadriven */
} /* namespace sgpp */
