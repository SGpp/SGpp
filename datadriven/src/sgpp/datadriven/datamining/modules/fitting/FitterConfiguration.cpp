// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "FitterConfiguration.hpp"

#include <string>

namespace sgpp {
namespace datadriven {

FitterConfiguration::FitterConfiguration()
    : gridConfig(),
      adaptivityConfig(),
      solverRefineConfig(),
      solverFinalConfig(),
      regularizationConfig(),
      lambda(0) {}

sgpp::datadriven::FitterConfiguration::~FitterConfiguration() {}

base::RegularGridConfiguration& FitterConfiguration::getGridConfig() { return gridConfig; }

base::AdpativityConfiguration& FitterConfiguration::getRefinementConfig() {
  return adaptivityConfig;
}

solver::SLESolverConfiguration& FitterConfiguration::getSolverRefineConfig() {
  return solverRefineConfig;
}

solver::SLESolverConfiguration& FitterConfiguration::getSolverFinalConfig() {
  return solverFinalConfig;
}

datadriven::RegularizationConfiguration& FitterConfiguration::getRegularizationConfig() {
  return regularizationConfig;
}

double FitterConfiguration::getLambda() { return lambda; }

void FitterConfiguration::setGridConfig(const base::RegularGridConfiguration& gridConfig) {
  this->gridConfig = gridConfig;
}

void FitterConfiguration::setRefinementConfig(
    const base::AdpativityConfiguration& adaptivityConfig) {
  this->adaptivityConfig = adaptivityConfig;
}

void FitterConfiguration::setSolverRefineConfig(
    const solver::SLESolverConfiguration& solverConfig) {
  this->solverRefineConfig = solverConfig;
}

void FitterConfiguration::setSolverFinalConfig(const solver::SLESolverConfiguration& solverConfig) {
  this->solverFinalConfig = solverConfig;
}

void FitterConfiguration::setRegularizationConfig(
    const datadriven::RegularizationConfiguration& regularizationConfig) {
  this->regularizationConfig = regularizationConfig;
}

void FitterConfiguration::setLambda(double lambda) { this->lambda = lambda; }

}  // namespace datadriven
}  // namespace sgpp
