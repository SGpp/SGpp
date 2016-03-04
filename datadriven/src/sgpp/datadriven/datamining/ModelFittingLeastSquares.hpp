// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include "ModelFittingBase.hpp"

#include "DataMiningConfigurationLeastSquares.hpp"

#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>
#include <sgpp/solver/SLESolver.hpp>

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
 private:
  base::OperationMultipleEval* kernel = nullptr;

  datadriven::OperationMultipleEvalConfiguration implementationConfiguration;

  DataMiningConfigurationLeastSquares configuration;

  std::shared_ptr<datadriven::DMSystemMatrixBase> systemMatrix;

  std::shared_ptr<solver::SLESolver> solver;

 protected:
  virtual datadriven::DMSystemMatrixBase* createSystemMatrix(base::DataMatrix& trainDataset,
                                                             double lambda);

 public:
  /**
   * Constructor
   *
   * @param config configuration
   */
  ModelFittingLeastSquares(sgpp::datadriven::DataMiningConfigurationLeastSquares config);

  /**
   * Destructor
   */
  virtual ~ModelFittingLeastSquares();

  // new grid and new dataset
  void fit(datadriven::Dataset& dataset) override;

  // reuse grid and assume old dataset
  // for grid refinement steps
  void refine() override;

  // reuse grid and new dataset
  // for online learning
  void update(datadriven::Dataset& dataset) override;

  void setImplementation(datadriven::OperationMultipleEvalConfiguration operationConfiguration) {
    this->implementationConfiguration = operationConfiguration;
  }
};
}
}
