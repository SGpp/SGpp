// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>

#include <sgpp/solver/SLESolver.hpp>

using sgpp::solver::SLESolver;
using sgpp::base::DataMatrix;

namespace sgpp {
namespace datadriven {

/**
 * This class implements standard sparse grid regression
 * with an Identity matrix as regularization operator.
 *
 * Furthermore this Learner provides support for several
 * vectorization approaches covering GPUs, CPUs and coprocessors.
 */
class ModelFittingLeastSquares : public ModelFittingBase {
 public:
  /**
   * Constructor
   *
   * @param config configuration
   */
  ModelFittingLeastSquares(const FitterConfigurationLeastSquares& config);

  /**
   * Destructor
   */
  virtual ~ModelFittingLeastSquares();

  // new grid and new dataset
  void fit(Dataset& dataset) override;

  // reuse grid and assume old dataset
  // for grid refinement steps
  void refine() override;

  // reuse grid and new dataset
  // for online learning
  void update(Dataset& dataset) override;
  //
  //  void setImplementation(datadriven::OperationMultipleEvalConfiguration operationConfiguration)
  //  {
  //    this->implementationConfiguration = operationConfiguration;
  //  }

 protected:
  virtual DMSystemMatrixBase* buildSystemMatrix(DataMatrix& trainDataset, double lambda);

  virtual SLESolver* buildSolver(FitterConfiguration& config);

  void configureSolver(FitterConfiguration& config, SLESolver& solver,
                       FittingSolverState solverState);

 private:
  FitterConfigurationLeastSquares config;
  std::shared_ptr<DMSystemMatrixBase> systemMatrix;
  std::shared_ptr<SLESolver> solver;
  OperationMultipleEvalConfiguration implementationConfig;
};
}
}
