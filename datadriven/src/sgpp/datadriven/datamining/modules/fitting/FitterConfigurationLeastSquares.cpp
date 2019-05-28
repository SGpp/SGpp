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
  FitterConfiguration::setupDefaults();
  // (Sebastian) Now all default-setup moved to parent class
}

void FitterConfigurationLeastSquares::readParams(const DataMiningConfigParser &parser) {
  setupDefaults();

  parser.getFitterGridConfig(gridConfig, gridConfig);
  // parser.getFitterAdaptivityConfig(adaptivityConfig, adaptivityConfig);
  parser.getFitterCoarseningConfig(coarseningConfig, coarseningConfig);
  parser.getFitterSolverRefineConfig(solverRefineConfig, solverRefineConfig);
  parser.getFitterSolverFinalConfig(solverFinalConfig, solverFinalConfig);
  parser.getFitterRegularizationConfig(regularizationConfig, regularizationConfig);
}
} /* namespace datadriven */
} /* namespace sgpp */
