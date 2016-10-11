/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FitterConfigurationLeastSquares.hpp
 *
 * Created on: Oct 10, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include "FitterConfiguration.hpp"

namespace sgpp {
namespace datadriven {

class FitterConfigurationLeastSquares : public FitterConfiguration {
 public:
  FitterConfigurationLeastSquares() : FitterConfiguration() { setupDefaults(); }

  FitterConfigurationLeastSquares(const DataMiningConfigParser& parser) : FitterConfiguration() {
    setupDefaults();
    parser.getFitterGridConfig(gridConfig, gridConfig);
    parser.getFitterAdaptivityConfig(adaptivityConfig, adaptivityConfig);
    parser.getFitterSolverRefineConfig(solverRefineConfig, solverRefineConfig);
    parser.getFitterSolverFinalConfig(solverFinalConfig, solverFinalConfig);
    parser.getFitterRegularizationConfig(regularizationConfig, regularizationConfig);
    parser.getFitterLambda(lambda, lambda);
  }

 private:
  void setupDefaults() {
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
    regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Laplace;

    lambda = 0.0;
  };
};

} /* namespace datadriven */
} /* namespace sgpp */
