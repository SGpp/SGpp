// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "DataMiningConfigurationLeastSquares.hpp"

#include <string>

#include "sgpp/base/tools/json/json_exception.hpp"

namespace sgpp {
namespace datadriven {

DataMiningConfigurationLeastSquares::DataMiningConfigurationLeastSquares()
    : DataMiningConfiguration() {
  // set default config
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
  regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Laplace;

  lambda = 0.0;
}

DataMiningConfigurationLeastSquares::DataMiningConfigurationLeastSquares(
    const std::string &fileName)
    : DataMiningConfiguration(fileName) {
  try {
    gridConfig.dim_ = 0;
    gridConfig.level_ = static_cast<int>((*this)["grid"]["level"].getUInt());
    gridConfig.type_ = stringToGridType((*this)["grid"]["type"].get());

    // configure adaptive refinement
    adaptivityConfig.maxLevelType_ = (*this)["refinement"]["maxLevelType"].getBool();
    adaptivityConfig.noPoints_ = (*this)["refinement"]["numPoints"].getUInt();
    adaptivityConfig.numRefinements_ = (*this)["refinement"]["numRefinements"].getUInt();
    adaptivityConfig.percent_ = (*this)["refinement"]["percent"].getDouble();
    adaptivityConfig.threshold_ = (*this)["refinement"]["threshold"].getDouble();

    // configure solver
    solverRefineConfig.type_ = stringToSolverType((*this)["solverRefine"]["type"].get());
    solverRefineConfig.maxIterations_ = (*this)["solverRefine"]["maxIterations"].getUInt();
    solverRefineConfig.eps_ = (*this)["solverRefine"]["eps"].getDouble();
    solverRefineConfig.threshold_ = (*this)["solverRefine"]["threshold"].getDouble();

    // configure solver
    solverFinalConfig.type_ = stringToSolverType((*this)["solverFinal"]["type"].get());
    solverFinalConfig.maxIterations_ = (*this)["solverFinal"]["maxIterations"].getUInt();
    solverFinalConfig.eps_ = (*this)["solverFinal"]["eps"].getDouble();
    solverFinalConfig.threshold_ = (*this)["solverFinal"]["threshold"].getDouble();

    // configure regularization
    regularizationConfig.regType_ =
        stringToRegularizationType((*this)["regularization"]["type"].get());

    lambda = (*this)["regularization"]["type"].getDouble();
  } catch (json::json_exception &e) {
    std::cout << e.what() << std::endl;
  }
}

base::RegularGridConfiguration DataMiningConfigurationLeastSquares::getGridConfig() {
  return gridConfig;
}

base::AdpativityConfiguration DataMiningConfigurationLeastSquares::getRefinementConfig() {
  return adaptivityConfig;
}

solver::SLESolverConfiguration DataMiningConfigurationLeastSquares::getSolverRefineConfig() {
  return solverRefineConfig;
}

solver::SLESolverConfiguration DataMiningConfigurationLeastSquares::getSolverFinalConfig() {
  return solverFinalConfig;
}

datadriven::RegularizationConfiguration
DataMiningConfigurationLeastSquares::getRegularizationConfig() {
  return regularizationConfig;
}

double DataMiningConfigurationLeastSquares::getLambda() { return lambda; }

void DataMiningConfigurationLeastSquares::setGridConfig(
    base::RegularGridConfiguration &gridConfig) {
  this->gridConfig = gridConfig;
}

void DataMiningConfigurationLeastSquares::setRefinementConfig(
    base::AdpativityConfiguration &adaptivityConfig) {
  this->adaptivityConfig = adaptivityConfig;
}

void DataMiningConfigurationLeastSquares::setSolverRefineConfig(
    solver::SLESolverConfiguration &solverConfig) {
  this->solverRefineConfig = solverConfig;
}

void DataMiningConfigurationLeastSquares::setSolverFinalConfig(
    solver::SLESolverConfiguration &solverConfig) {
  this->solverFinalConfig = solverConfig;
}

void DataMiningConfigurationLeastSquares::setRegularizationConfig(
    datadriven::RegularizationConfiguration &regularizationConfig) {
  this->regularizationConfig = regularizationConfig;
}

void DataMiningConfigurationLeastSquares::setLambda(double lambda) { this->lambda = lambda; }

DataMiningConfiguration *DataMiningConfigurationLeastSquares::clone() {
  DataMiningConfiguration *clone = new DataMiningConfigurationLeastSquares(*this);
  return clone;
}

}  // namespace datadriven
}  // namespace sgpp
