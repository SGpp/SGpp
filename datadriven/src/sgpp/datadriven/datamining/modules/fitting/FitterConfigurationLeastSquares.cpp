/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FitterConfigurationLeastSquares.cpp
 *
 *  Created on: 25.01.2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp>

namespace sgpp {
namespace datadriven {

FitterConfiguration *FitterConfigurationLeastSquares::clone() const {
  return new FitterConfigurationLeastSquares(*this);
}

void FitterConfigurationLeastSquares::setupDefaults() {
  // configure initial grid
  gridConfig.dim_ = 0;
  gridConfig.level_ = 2;
  gridConfig.type_ = sgpp::base::GridType::Linear;
  gridConfig.maxDegree_ = 1;
  gridConfig.boundaryLevel_ = 0;

  // configure adaptive refinement
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.noPoints_ = 0;
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.percent_ = 100.0;
  adaptivityConfig.threshold_ = 0.0,

      // configure solver
      solverRefineConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverRefineConfig.maxIterations_ = 100;
  solverRefineConfig.eps_ = 1e-10;
  solverRefineConfig.threshold_ = 1e-10;

  // configure solver
  solverFinalConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverFinalConfig.maxIterations_ = 100;
  solverFinalConfig.eps_ = 1e-10;
  solverFinalConfig.threshold_ = 1e-10;

  // configure regularization
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 0.0;
  regularizationConfig.l1Ratio_ = 0.0;
  regularizationConfig.exponentBase_ = 0.0;
}

void FitterConfigurationLeastSquares::readParams(const DataMiningConfigParser &parser) {
  setupDefaults();

  parser.getFitterGridConfig(gridConfig, gridConfig);
  //parser.getFitterAdaptivityConfig(adaptivityConfig, adaptivityConfig);
  parser.getFitterCoarseningConfig(coarseningConfig, coarseningConfig);
  parser.getFitterSolverRefineConfig(solverRefineConfig, solverRefineConfig);
  parser.getFitterSolverFinalConfig(solverFinalConfig, solverFinalConfig);
  parser.getFitterRegularizationConfig(regularizationConfig, regularizationConfig);
}
} /* namespace datadriven */
} /* namespace sgpp */
