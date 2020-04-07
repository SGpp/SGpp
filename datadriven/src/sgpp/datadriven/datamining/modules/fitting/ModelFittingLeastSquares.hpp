// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>
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
class ModelFittingLeastSquares : public ModelFittingBaseSingleGrid {
 public:
  /**
   * Constructor
   *
   * @param config configuration object that specifies grid, refinement, and regularization
   */
  explicit ModelFittingLeastSquares(const FitterConfigurationLeastSquares &config);

  /**
   * Fit the grid to the given dataset by determining the weights of the initial grid by a least
   * squares approach.
   * @param dataset the training dataset that is used to fit the model.
   */
  void fit(Dataset &dataset) override;
  void fit(Dataset &, Dataset &) override {
    throw base::application_exception("This model requires a single input dataset");
  }

  /**
   * Improve accuracy of the fit on the given training data by adaptive refinement of the grid and
   * recalculate weights.
   * @return true if refinement could be performed based on the refinement configuration, else
   * false.
   */
  bool adapt() override;

  void update(Dataset &dataset) override;
  void update(Dataset &, Dataset &) override {
    throw base::application_exception("This model requires a single input dataset");
  }

  /**
   * Evaluate the fitted regression model at a single data point - requires a trained grid.
   * @param sample vector with the coordinates in all dimensions of that sample.
   * @return evaluation of the trained grid.
   */
  double evaluate(const DataVector &sample) override;

  /**
   * Evaluate the fitted model on a set of data points - requires a trained grid.
   * @param samples matrix where each row represents a sample and the columns contain the
   * coordinates in all dimensions of that sample.
   * @param results vector where each row will contain the evaluation of the respective sample on
   * the current model.
   */
  void evaluate(DataMatrix &samples, DataVector &results) override;

  /**
   * Resets the state of the entire model
   */
  void reset() override;

  /**
   * Resets any trained representations of the model, but does not reset the entire state.
   */
  void resetTraining() override {
    throw sgpp::base::not_implemented_exception(
        "ModelFittingLeastSquares::resetTraining() is not implemented!");
  }

  /**
   * Should compute some kind of Residual to evaluate the fit of the model.
   *
   * In the case of density estimation, this is
   * || R * alpha_lambda - b_val ||_2
   *
   * This is useful for unsupervised learning models, where normal evaluation cannot be used as
   * there are no targets.
   *
   * @param validationData Matrix for validation data
   *
   * @returns the residual score
   */
  double computeResidual(DataMatrix &validationData) const override {
    throw sgpp::base::not_implemented_exception(
        "ModelFittingLeastSquares::computeResidual() is not implemented!");
  }

  /**
   * Updates the regularization parameter lambda of the underlying model.
   *
   * @param lambda the new lambda parameter
   */
  void updateRegularization(double lambda) override {
    throw sgpp::base::not_implemented_exception(
        "ModelFittingLeastSquares::updateRegularization() is not implemented!");
  }

 private:
  /**
   * Count the amount of refinement operations performed on the current dataset.
   */
  size_t refinementsPerformed;

  // TODO(lettrich): grid and train dataset as well as OperationMultipleEvalConfiguration should be
  // const.
  /**
   * Factory function to build the System matrix for least squares regression with identity as
   * regularization.
   */
  DMSystemMatrixBase *buildSystemMatrix(Grid &grid, DataMatrix &trainDataset, double lambda,
                                        OperationMultipleEvalConfiguration &config) const;

  /**
   * based on the current dataset and grid, assemble a system of linear equations and solve for the
   * hierarchical surplus vector alpha.
   * @param solverConfig: Configuration of the SLESolver (refinement, or final solver).
   * @param alpha: Reference to a data vector where hierarchical surpluses will be stored into. Make
   * sure the vector size is equal to the amount of grid points.
   */
  void assembleSystemAndSolve(const SLESolverConfiguration &solverConfig, DataVector &alpha) const;
};
} /* namespace datadriven */
} /* namespace sgpp */
