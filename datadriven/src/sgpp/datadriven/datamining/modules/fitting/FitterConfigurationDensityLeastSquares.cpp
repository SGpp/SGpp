// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityLeastSquares.hpp>

namespace sgpp {
namespace datadriven {

FitterConfiguration *FitterConfigurationDensityLeastSquares::clone() const {
  return new FitterConfigurationDensityLeastSquares(*this);
}

void FitterConfigurationDensityLeastSquares::setupDefaults() {
  FitterConfiguration::setupDefaults();
  // (Sebastian) Now all default-setup moved to parent class
}

void FitterConfigurationDensityLeastSquares::readParams(const DataMiningConfigParser &parser) {
  setupDefaults();

  parser.getFitterGridConfig(gridConfig, gridConfig);
  parser.getFitterAdaptivityConfig(adaptivityConfig, adaptivityConfig);
  parser.getFitterSolverRefineConfig(solverRefineConfig, solverRefineConfig);
  parser.getFitterSolverFinalConfig(solverFinalConfig, solverFinalConfig);
  parser.getFitterRegularizationConfig(regularizationConfig, regularizationConfig);
  // add the density estimation config
  parser.getFitterDensityEstimationConfig(densityEstimationConfig, densityEstimationConfig);
}
} /* namespace datadriven */
} /* namespace sgpp */
