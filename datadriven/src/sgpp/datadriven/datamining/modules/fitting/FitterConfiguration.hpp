// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/solver/TypesSolver.hpp>

namespace sgpp {
namespace datadriven {

enum class FitterType { RegressionLeastSquares };

class FitterConfiguration {
 public:
  FitterConfiguration();
  FitterConfiguration(const FitterConfiguration& rhs) = default;
  FitterConfiguration(FitterConfiguration&& rhs) = default;
  FitterConfiguration& operator=(const FitterConfiguration& rhs) = default;
  FitterConfiguration& operator=(FitterConfiguration&& rhs) = default;
  virtual ~FitterConfiguration() = default;

  const base::RegularGridConfiguration& getGridConfig() const;
  const base::AdpativityConfiguration& getRefinementConfig() const;
  const solver::SLESolverConfiguration& getSolverRefineConfig() const;
  const solver::SLESolverConfiguration& getSolverFinalConfig() const;
  const datadriven::RegularizationConfiguration& getRegularizationConfig() const;
  const datadriven::OperationMultipleEvalConfiguration& getMultipleEvalConfig() const;

  base::RegularGridConfiguration& getGridConfig();
  base::AdpativityConfiguration& getRefinementConfig();
  solver::SLESolverConfiguration& getSolverRefineConfig();
  solver::SLESolverConfiguration& getSolverFinalConfig();
  datadriven::RegularizationConfiguration& getRegularizationConfig();
  datadriven::OperationMultipleEvalConfiguration& getMultipleEvalConfig();

  double getLambda();
  void setLambda(double lambda);

 protected:
  virtual void setupDefaults() = 0;
  base::RegularGridConfiguration gridConfig;
  base::AdpativityConfiguration adaptivityConfig;
  solver::SLESolverConfiguration solverRefineConfig;
  solver::SLESolverConfiguration solverFinalConfig;
  datadriven::RegularizationConfiguration regularizationConfig;
  double lambda;
  datadriven::OperationMultipleEvalConfiguration multipleEvalConfig;
};

}  // namespace base
}  // namespace sgpp
