// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>

namespace sgpp {
namespace datadriven {

class FitterConfiguration {
 public:
  FitterConfiguration();

  // we want to make this an abstract class.
  virtual ~FitterConfiguration() = 0;

  base::RegularGridConfiguration& getGridConfig();
  base::AdpativityConfiguration& getRefinementConfig();
  solver::SLESolverConfiguration& getSolverRefineConfig();
  solver::SLESolverConfiguration& getSolverFinalConfig();
  datadriven::RegularizationConfiguration& getRegularizationConfig();
  double getLambda();

  void setGridConfig(const base::RegularGridConfiguration& gridConfig);
  void setRefinementConfig(const base::AdpativityConfiguration& adaptivityConfig);
  void setSolverRefineConfig(const solver::SLESolverConfiguration& solverRefineConfig);
  void setSolverFinalConfig(const solver::SLESolverConfiguration& solverFinalConfig);
  void setRegularizationConfig(const datadriven::RegularizationConfiguration& regularizationConfig);
  void setLambda(double lambda);

 protected:
  base::RegularGridConfiguration gridConfig;
  base::AdpativityConfiguration adaptivityConfig;
  solver::SLESolverConfiguration solverRefineConfig;
  solver::SLESolverConfiguration solverFinalConfig;
  datadriven::RegularizationConfiguration regularizationConfig;
  double lambda;
};

}  // namespace base
}  // namespace sgpp
