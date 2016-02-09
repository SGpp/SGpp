// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "DataMiningConfigurationLeastSquares.hpp"

#include <string>

namespace SGPP {
namespace datadriven {

DataMiningConfigurationLeastSquares::DataMiningConfigurationLeastSquares() :
		DataMiningConfiguration() {
	// set default config
	gridConfig.dim_ = 0;
	gridConfig.level_ = 2;
	gridConfig.type_ = SGPP::base::GridType::Linear;
	gridConfig.maxDegree = 1;
	gridConfig.boundaryLevel = 0;

	// configure adaptive refinement
	adaptivityConfig.maxLevelType_ = false;
	adaptivityConfig.noPoints_ = 0;
	adaptivityConfig.numRefinements_ = 0;
	adaptivityConfig.percent_ = 100.0;
	adaptivityConfig.threshold_ = 0.0,

	// configure solver
	solverRefineConfig.type_ = SGPP::solver::SLESolverType::CG;
	solverRefineConfig.maxIterations_ = 100;
	solverRefineConfig.eps_ = 1e-10;
	solverRefineConfig.threshold_ = 1e-10;

	// configure solver
	solverFinalConfig.type_ = SGPP::solver::SLESolverType::CG;
	solverFinalConfig.maxIterations_ = 100;
	solverFinalConfig.eps_ = 1e-10;
	solverFinalConfig.threshold_ = 1e-10;

	// configure regularization
	regularizationConfig.regType_ =
			SGPP::datadriven::RegularizationType::Laplace;

	lambda = 0.0;
}


DataMiningConfigurationLeastSquares::DataMiningConfigurationLeastSquares(
		const std::string& fileName) :
		DataMiningConfiguration(fileName) {
}

base::RegularGridConfiguration &DataMiningConfigurationLeastSquares::getGridConfig() {
	return gridConfig;
}

base::AdpativityConfiguration &DataMiningConfigurationLeastSquares::getRefinementConfig() {
	return adaptivityConfig;
}

solver::SLESolverConfiguration &DataMiningConfigurationLeastSquares::getSolverRefineConfig() {
	return solverRefineConfig;
}

solver::SLESolverConfiguration &DataMiningConfigurationLeastSquares::getSolverFinalConfig() {
	return solverFinalConfig;
}

datadriven::RegularizationConfiguration &DataMiningConfigurationLeastSquares::getRegularizationConfig() {
	return regularizationConfig;
}

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

DataMiningConfiguration* DataMiningConfigurationLeastSquares::clone() {
	DataMiningConfiguration* clone = new DataMiningConfigurationLeastSquares(*this);
	return clone;
}

}// namespace base
}  // namespace SGPP

