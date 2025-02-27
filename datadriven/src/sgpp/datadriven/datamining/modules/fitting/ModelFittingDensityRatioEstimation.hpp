// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrixTwoDatasets.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/solver/SLESolver.hpp>

using sgpp::solver::SLESolver;
using sgpp::base::DataMatrix;
using sgpp::base::Grid;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven {

/**
 * Fitter object that encapsulates the usage of sparse grid based density ratio estimation, based
 * on the regression model with identity as regularization.
 *
 * Allows usage of different grids, different solvers and different regularization techniques
 * based on the provided configuration objects.
 */
class ModelFittingDensityRatioEstimation : public ModelFittingBaseSingleGrid {
 public:
  /**
   * Constructor
   *
   * @param config configuration object that specifies grid, refinement, and regularization
   */
  explicit ModelFittingDensityRatioEstimation(const FitterConfigurationLeastSquares &config);

  /**
   * Fit the grid to the given dataset by determining the weights of the initial grid by a least
   * squares approach.
   * @param newDatasetP the training dataset of first density that is used to fit the model
   * @param newDatasetQ the training dataset of first density that is used to fit the model
   */
  void fit(Dataset &newDatasetP, Dataset &newDatasetQ) override;
  void fit(Dataset &) override {
    throw base::application_exception("This model requires two input datasets");
  }

  /**
   * Improve accuracy of the fit on the given training data by adaptive refinement of the grid and
   * recalculate weights.
   * @return true if refinement could be performed based on the refinement configuration, else
   * false.
   */
  bool adapt() override;

  void update(Dataset &datasetP, Dataset &datasetQ) override;
  void update(Dataset &) override {
    throw base::application_exception("This model requires two input datasets");
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
   * Computes an approximation of the least squares loss between true and estimated density ratio.
   * Implemented as: 1 / M_q * Sum f(x_q)^2 - 2 / M_p * Sum f(x_p)
   * @param samplesP samples of first dataset to evaluate against
   * @param samplesQ samples of second dataset to evaluate against
   * @return LS loss approximation
   */
  double LeastSquaresLossApprox(DataMatrix &samplesP, DataMatrix &samplesQ) override;

  /**
   * Computes an approximation of the Kullback-Leibler divergence using the density ratio.
   * Implemented as: 1 / M_p * Sum log f(x_p)
   * @param samplesP samples of first dataset to evaluate against
   * @return PE divergence approximation
   */
  double KLDivergenceApprox(DataMatrix &samplesP) override;

  /**
   * Computes an approximation of the Pearson divergence using the density ratio.
   * Implemented as: -1/2 * 1 / M_q * Sum f(x_q)^2 + 1 / M_p * Sum f(x_p) - 1/2
   * @param samplesP samples of first dataset to evaluate against
   * @param samplesQ samples of second dataset to evaluate against
   * @return PE divergence approximation
   */
  double PEDivergenceApprox(DataMatrix &samplesP, DataMatrix &samplesQ) override;

  /**
   * Resets the state of the entire model
   */
  void reset() override;

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
        "ModelFittingDensityRatioEstimation::computeResidual() is not implemented!");
  }

  /**
   * Updates the regularization parameter lambda of the underlying model.
   *
   * @param lambda the new lambda parameter
   */
  void updateRegularization(double lambda) override {
    throw sgpp::base::not_implemented_exception(
        "ModelFittingDensityRatioEstimation::updateRegularization() is not implemented!");
  }

  /**
   * Resets any trained representations of the model, but does not reset the entire state.
   */
  void resetTraining() override {
    throw sgpp::base::not_implemented_exception(
        "ModelFittingDensityRatioEstimation::resetTraining() is not implemented!");
  }

 private:
  /**
   * Returns the refinement functor suitable for the model settings.
   * @return pointer to a refinement functor that suits the model settings
   */
  std::unique_ptr<sgpp::base::RefinementFunctor> getRefinementFunctor();

  /**
   * Returns the refinement functor suitable for the model settings.
   * @return pointer to a coarsening functor that suits the model settings
   */
  std::unique_ptr<sgpp::base::CoarseningFunctor> getCoarseningFunctor();

  /**
   * Count the amount of refinement operations performed on the current dataset.
   */
  size_t refinementsPerformed;

  /**
   * Initial number of grid points
   */
  size_t initialGridSize;

  // TODO(lettrich): grid and train dataset as well as OperationMultipleEvalConfiguration should
  // be const.
  /**
   * Factory function to build the System matrix for density ratio estimation with identity as
   * regularization
   */
  DMSystemMatrixTwoDatasets *buildSystemMatrix(Grid &grid, DataMatrix &trainDatasetP,
                                               DataMatrix &trainDatasetQ, double lambda,
                                               OperationMultipleEvalConfiguration &config) const;

  /**
   * Based on the current dataset and grid, assemble a system of linear equations and solve for
   * the hierarchical surplus vector alpha.
   * @param solverConfig: Configuration of the SLESolver (refinement, or final solver).
   * @param alpha: Reference to a data vector where hierarchical surpluses will be stored into.
   * Make sure the vector size is equal to the amount of grid points.
   */
  void assembleSystemAndSolve(const SLESolverConfiguration &solverConfig, DataVector &alpha) const;
};
} /* namespace datadriven */
} /* namespace sgpp */
