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
using sgpp::base::Grid;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven {

// TODO(lettrich): allow different refinement techniques.
/**
 * Fitter object that encapsulates the usage of sparse grid based regression with identity as
 * regularization.
 *
 * Allows usage of different grids, different solvers and different regularization techniques based
 * on the provided configuration objects.
 */
class ModelFittingLeastSquares : public ModelFittingBase {
 public:
  /**
   * Constructor
   *
   * @param config configuration object that specifies grid, refinement, and regularization
   */
  explicit ModelFittingLeastSquares(const FitterConfigurationLeastSquares& config);

  /**
   * Fit the grid to the given dataset by determining the weights of the initial grid by a least
   * squares approach.
   * @param dataset the training dataset that is used to fit the model.
   */
  void fit(Dataset& dataset) override;

  /**
   * Improve accuracy of the fit on the given training data by adaptive refinement of the grid and
   * recalculate weights.
   */
  void refine() override;

  void update(Dataset& dataset) override;

  /**
   * Evaluate the fitted regression model at a single data point - requires a trained grid.
   * @param sample vector with the coordinates in all dimensions of that sample.
   * @return evaluation of the trained grid.
   */
  double evaluate(const DataVector& sample) const override;

  /**
   * Evaluate the fitted model on a set of data points - requires a trained grid.
   * @param samples matrix where each row represents a sample and the columns contain the
   * coordinates in all dimensions of that sample.
   * @param results vector where each row will contain the evaluation of the respective sample on
   * the current model.
   */
  void evaluate(DataMatrix& samples, DataVector& results) override;

 private:
  // TODO(lettrich): grid and train dataset as well as OperationMultipleEvalConfiguration should be
  // const.
  /**
   * Factory function to build the System matrix for least squares regression with identity as
   * regularization.
   */
  DMSystemMatrixBase* buildSystemMatrix(Grid& grid, DataMatrix& trainDataset, double lambda,
                                        OperationMultipleEvalConfiguration& config) const;
};
} /* namespace datadriven */
} /* namespace sgpp */
