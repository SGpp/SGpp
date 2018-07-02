/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FitterConfigurationClassification.cpp
 *
 *  Created on: Jul 2, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationClassification.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>

namespace sgpp {
namespace datadriven {

FitterConfiguration* FitterConfigurationClassification::clone() const {
  return new FitterConfigurationClassification(*this);
}

void FitterConfigurationClassification::setupDefaults() {
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
}

void FitterConfigurationClassification::readParams(const DataMiningConfigParser& parser) {
  setupDefaults();

  parser.getFitterGridConfig(gridConfig, gridConfig);
  parser.getFitterAdaptivityConfig(adaptivityConfig, adaptivityConfig);
  parser.getFitterSolverRefineConfig(solverRefineConfig, solverRefineConfig);
  parser.getFitterSolverFinalConfig(solverFinalConfig, solverFinalConfig);
  parser.getFitterRegularizationConfig(regularizationConfig, regularizationConfig);
  parser.getFitterDatabaseConfig(databaseConfig, databaseConfig);
  parser.getFitterLearnerConfig(learnerConfig, learnerConfig);
}
} /* namespace datadriven */
} /* namespace sgpp */
