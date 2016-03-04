// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include "DataMiningConfiguration.hpp"
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>

#include <vector>
#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

class DataMiningConfigurationLeastSquares : public DataMiningConfiguration {
  friend class ModelFittingDensityEstimation;

 private:
  base::RegularGridConfiguration gridConfig;
  base::AdpativityConfiguration adaptivityConfig;
  solver::SLESolverConfiguration solverRefineConfig;
  solver::SLESolverConfiguration solverFinalConfig;
  datadriven::RegularizationConfiguration regularizationConfig;
  double lambda;

 public:
  DataMiningConfigurationLeastSquares();

  explicit DataMiningConfigurationLeastSquares(const std::string &fileName);

  base::RegularGridConfiguration getGridConfig();
  base::AdpativityConfiguration getRefinementConfig();
  solver::SLESolverConfiguration getSolverRefineConfig();
  solver::SLESolverConfiguration getSolverFinalConfig();
  datadriven::RegularizationConfiguration getRegularizationConfig();
  double getLambda();

  void setGridConfig(base::RegularGridConfiguration &gridConfig);
  void setRefinementConfig(base::AdpativityConfiguration &adaptivityConfig);
  void setSolverRefineConfig(solver::SLESolverConfiguration &solverRefineConfig);
  void setSolverFinalConfig(solver::SLESolverConfiguration &solverFinalConfig);
  void setRegularizationConfig(datadriven::RegularizationConfiguration &regularizationConfig);
  void setLambda(double lambda);

  virtual DataMiningConfiguration *clone() override;
};

}  // namespace base
}  // namespace sgpp
