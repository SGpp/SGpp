// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/DataMiningConfigurationDensityEstimation.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <string>

namespace SGPP {
namespace datadriven {

DataMiningConfigurationDensityEstimation::DataMiningConfigurationDensityEstimation(): SGPP::datadriven::DataMiningConfiguration() {
  // set default config
  gridConfig.dim_ = 1; // dataset.getDimension();
  gridConfig.level_ = 2;
  gridConfig.type_ = SGPP::base::GridType::Linear;

  // configure adaptive refinement
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.noPoints_ = 0;

  // configure solver
  solverConfig.type_ = SGPP::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 100;
  solverConfig.eps_ = 1e-10;
  solverConfig.threshold_ = 1e-10;

  // configure regularization
  regularizationConfig.regType_ = SGPP::datadriven::RegularizationType::Laplace;

  // configure learner
  learnerSGDEConfig.doCrossValidation_ = false;
  learnerSGDEConfig.kfold_ = 10;
  learnerSGDEConfig.lambdaStart_ = 1e-1;
  learnerSGDEConfig.lambdaEnd_ = 1e-10;
  learnerSGDEConfig.lambdaSteps_ = 10;
  learnerSGDEConfig.logScale_ = true;
  learnerSGDEConfig.shuffle_ = true;
  learnerSGDEConfig.seed_ = 1234567;
  learnerSGDEConfig.silent_ = true;
}

DataMiningConfigurationDensityEstimation::DataMiningConfigurationDensityEstimation(const std::string& fileName):
    SGPP::datadriven::DataMiningConfiguration(fileName) {
  // initialize structs from file
  // configure grid
  try {
    gridConfig.dim_ = 1; // dataset.getDimension();
    gridConfig.level_ = (*this)["grid"]["level"].getUInt();
    gridConfig.type_ = stringToGridType((*this)["grid"]["type"].get());

    // configure adaptive refinement
    adaptivityConfig.numRefinements_ = (*this)["refinement"]["numSteps"].getUInt();
    adaptivityConfig.noPoints_ = (*this)["refinement"]["numPoints"].getUInt();

    // configure solver
    solverConfig.type_ = stringToSolverType((*this)["solver"]["type"].get());
    solverConfig.maxIterations_ = (*this)["solver"]["maxIterations"].getUInt();
    solverConfig.eps_ = (*this)["solver"]["eps"].getDouble();
    solverConfig.threshold_ = (*this)["solver"]["threshold"].getDouble();

    // configure regularization
    regularizationConfig.regType_ = stringToRegularizationType((*this)["regularization"]["type"].get());

    // configure learner
    learnerSGDEConfig.doCrossValidation_ = (*this)["crossValidation"]["doCrossValidation"].getBool();
    learnerSGDEConfig.kfold_ = (*this)["crossValidation"]["kfold"].getUInt();
    learnerSGDEConfig.lambdaStart_ = (*this)["crossValidation"]["lambdaStart"].getDouble();
    learnerSGDEConfig.lambdaEnd_ = (*this)["crossValidation"]["lambdaEnd"].getDouble();
    learnerSGDEConfig.lambdaSteps_ = (*this)["crossValidation"]["lambdaSteps"].getUInt();
    learnerSGDEConfig.logScale_ = (*this)["crossValidation"]["logScale"].getBool();
    learnerSGDEConfig.shuffle_ = (*this)["crossValidation"]["shuffle"].getBool();
    learnerSGDEConfig.seed_ = (*this)["crossValidation"]["seed"].getUInt();
    learnerSGDEConfig.silent_ = (*this)["crossValidation"]["verbose"].getBool();
  } catch (SGPP::base::application_exception& e) {
    std::cout << e.what() << std::endl;
  }
}

DataMiningConfiguration* DataMiningConfigurationDensityEstimation::clone() {
	DataMiningConfiguration* clone = new DataMiningConfigurationDensityEstimation(*this);
  return clone;
}

}  // namespace base
}  // namespace SGPP

