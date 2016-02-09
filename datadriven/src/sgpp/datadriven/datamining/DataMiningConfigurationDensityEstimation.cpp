// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/DataMiningConfigurationDensityEstimation.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <string>

namespace SGPP {
namespace datadriven {

DataMiningConfigurationDensityEstimation::DataMiningConfigurationDensityEstimation() :
		SGPP::datadriven::DataMiningConfiguration() {
	// set default config
	gridConfig.dim_ = 0;
	gridConfig.level_ = 2;
	gridConfig.type_ = SGPP::base::GridType::Linear;
	gridConfig.maxDegree_ = 1;
	gridConfig.boundaryLevel_ = 0;

	// configure adaptive refinement
	adaptivityConfig.numRefinements_ = 0;
	adaptivityConfig.noPoints_ = 0;

	// configure solver
	solverConfig.type_ = SGPP::solver::SLESolverType::CG;
	solverConfig.maxIterations_ = 100;
	solverConfig.eps_ = 1e-10;
	solverConfig.threshold_ = 1e-10;

	// configure regularization
	regularizationConfig.regType_ =
			SGPP::datadriven::RegularizationType::Laplace;

	// configure learner
	sgdeConfig.doCrossValidation_ = false;
	sgdeConfig.kfold_ = 10;
	sgdeConfig.lambdaStart_ = 1e-1;
	sgdeConfig.lambdaEnd_ = 1e-10;
	sgdeConfig.lambdaSteps_ = 10;
	sgdeConfig.logScale_ = true;
	sgdeConfig.shuffle_ = true;
	sgdeConfig.seed_ = 1234567;
	sgdeConfig.silent_ = true;
}

DataMiningConfigurationDensityEstimation::DataMiningConfigurationDensityEstimation(
		const std::string& fileName) :
		SGPP::datadriven::DataMiningConfiguration(fileName) {
	// initialize structs from file
	// configure grid
	try {
		gridConfig.dim_ = 0;
		gridConfig.level_ = (*this)["grid"]["level"].getUInt();
		gridConfig.type_ = stringToGridType((*this)["grid"]["type"].get());

		// configure adaptive refinement
		adaptivityConfig.numRefinements_ =
				(*this)["refinement"]["numSteps"].getUInt();
		adaptivityConfig.noPoints_ =
				(*this)["refinement"]["numPoints"].getUInt();

		// configure solver
		solverConfig.type_ = stringToSolverType(
				(*this)["solver"]["type"].get());
		solverConfig.maxIterations_ =
				(*this)["solver"]["maxIterations"].getUInt();
		solverConfig.eps_ = (*this)["solver"]["eps"].getDouble();
		solverConfig.threshold_ = (*this)["solver"]["threshold"].getDouble();

		// configure regularization
		regularizationConfig.regType_ = stringToRegularizationType(
				(*this)["regularization"]["type"].get());

		// configure learner
		sgdeConfig.doCrossValidation_ =
				(*this)["crossValidation"]["doCrossValidation"].getBool();
		sgdeConfig.kfold_ = (*this)["crossValidation"]["kfold"].getUInt();
		sgdeConfig.lambdaStart_ =
				(*this)["crossValidation"]["lambdaStart"].getDouble();
		sgdeConfig.lambdaEnd_ =
				(*this)["crossValidation"]["lambdaEnd"].getDouble();
		sgdeConfig.lambdaSteps_ =
				(*this)["crossValidation"]["lambdaSteps"].getUInt();
		sgdeConfig.logScale_ =
				(*this)["crossValidation"]["logScale"].getBool();
		sgdeConfig.shuffle_ = (*this)["crossValidation"]["shuffle"].getBool();
		sgdeConfig.seed_ = (*this)["crossValidation"]["seed"].getUInt();
		sgdeConfig.silent_ = (*this)["crossValidation"]["verbose"].getBool();
	} catch (SGPP::base::application_exception& e) {
		std::cout << e.what() << std::endl;
	}
}

base::RegularGridConfiguration &DataMiningConfigurationDensityEstimation::getGridConfig() {
	return gridConfig;
}

base::AdpativityConfiguration &DataMiningConfigurationDensityEstimation::getRefinementConfig() {
	return adaptivityConfig;
}

solver::SLESolverConfiguration &DataMiningConfigurationDensityEstimation::getSolverConfig() {
	return solverConfig;
}

datadriven::RegularizationConfiguration &DataMiningConfigurationDensityEstimation::getRegularizationConfig() {
	return regularizationConfig;
}

datadriven::DataMiningConfigurationDensityEstimationType &DataMiningConfigurationDensityEstimation::getSGDEConfig() {
	return sgdeConfig;
}

void DataMiningConfigurationDensityEstimation::setGridConfig(
		base::RegularGridConfiguration &gridConfig) {
	this->gridConfig = gridConfig;
}

void DataMiningConfigurationDensityEstimation::setRefinementConfig(
		base::AdpativityConfiguration &adaptivityConfig) {
	this->adaptivityConfig = adaptivityConfig;
}

void DataMiningConfigurationDensityEstimation::setSolverConfig(
		solver::SLESolverConfiguration &solverConfig) {
	this->solverConfig = solverConfig;
}

void DataMiningConfigurationDensityEstimation::setRegularizationConfig(
		datadriven::RegularizationConfiguration &regularizationConfig) {
	this->regularizationConfig = regularizationConfig;
}

void DataMiningConfigurationDensityEstimation::setSGDEConfig(
		datadriven::DataMiningConfigurationDensityEstimationType &sgdeConfig) {
	this->sgdeConfig = sgdeConfig;
}

DataMiningConfiguration* DataMiningConfigurationDensityEstimation::clone() {
	DataMiningConfiguration* clone =
			new DataMiningConfigurationDensityEstimation(*this);
	return clone;
}

}  // namespace base
}  // namespace SGPP

