// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>

#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::Grid;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven {

// TODO(lettrich): allow different refinement techniques.
/**
 * Fitter object that encapsulates the usage of sparse grid density difference estimation with
 * identity as regularization.
 *
 * Allows usage of different grids, different solvers and different regularization techniques based
 * on the provided configuration objects.
 */
class ModelFittingDensityDifferenceEstimationCG : public ModelFittingDensityEstimation {
 public:
  /**
   * Constructor
   *
   * @param config configuration object that specifies grid, refinement, and regularization
   */
  explicit ModelFittingDensityDifferenceEstimationCG(
      const FitterConfigurationDensityEstimation& config);

  /**
   * Fit the grid to the given dataset by determining the weights of the initial grid by the SGDE
   * approach.
   * @param datasetP the first training dataset that is used to fit the model.
   * @param datasetQ the second training dataset that is used to fit the model.
   */
  void fit(Dataset& datasetP, Dataset& datasetQ) override;
  void fit(Dataset& dataset) override {
    throw base::application_exception("This model requires two input datasets");
  }

  /**
   * Fit the grid to the given dataset by determining the weights of the initial grid by the SGDE
   * approach. Requires only data samples and no targets (since those are irrelevant for the density
   * estimation whatsoever)
   * @param datasetP the first training dataset that is used to fit the model.
   * @param datasetQ the second training dataset that is used to fit the model.
   */
  void fit(DataMatrix& datasetP, DataMatrix& datasetQ) override;
  void fit(DataMatrix& dataset) override {
    throw base::application_exception("This model requires two input datasets");
  }

  /**
   * Performs a refinement given the new grid size and the points to coarsened
   * @param newNoPoints the grid size after refinement and coarsening
   * @param deletedGridPoints a list of indexes for grid points that will be removed
   * @return if the grid was refined (true)
   */
  bool adapt(size_t newNoPoints, std::vector<size_t>& deletedGridPoints) override;

  void update(Dataset& datasetP, Dataset& datasetQ) override;
  void update(Dataset& dataset) override {
    throw base::application_exception("This model requires two input datasets");
  }

  /**
   * Updates the model based on new data samples (streaming, batch learning). Requires only the data
   * samples and no targets (since those are irrelevant for the density estimation whatsoever)
   * @param samplesP the new data samples for first dataset
   * @param samplesQ the new data samples for second dataset
   */
  void update(DataMatrix& samplesP, DataMatrix& samplesQ) override;
  void update(DataMatrix& samples) override {
    throw base::application_exception("This model requires two input datasets");
  }

  /**
   * Evaluate the fitted density at a single data point - requires a trained grid.
   * @param sample vector with the coordinates in all dimensions of that sample.
   * @return evaluation of the trained grid.
   */
  double evaluate(const DataVector& sample) override;

  /**
   * Evaluate the fitted density on a set of data points - requires a trained grid.
   * @param samples matrix where each row represents a sample and the columns contain the
   * coordinates in all dimensions of that sample.
   * @param results vector where each row will contain the evaluation of the respective sample on
   * the current model.
   */
  void evaluate(DataMatrix& samples, DataVector& results) override;

  /**
   * Resets the state of the entire model
   */
  void reset() override;

  /**
   * Resets any trained representations of the model, but does not reset the entire state.
   */
  void resetTraining() override {
    throw sgpp::base::not_implemented_exception(
        "ModelFittingDensityEstimationCG::resetTraining() is not implemented!");
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
  double computeResidual(DataMatrix& validationData) const override {
    throw sgpp::base::not_implemented_exception(
        "ModelFittingDensityEstimationCG::computeResidual() is not "
        "implemented!");
  }

  /**
   * Updates the regularization parameter lambda of the underlying model.
   *
   * @param lambda the new lambda parameter
   */
  void updateRegularization(double lambda) override {
    throw sgpp::base::not_implemented_exception(
        "ModelFittingDensityEstimationCG::updateRegularization() is not "
        "implemented!");
  }

 private:
  /**
   * Creates the regularization operation matrix for the model settings.
   * @param grid the underlying sparse grid
   * @return the regularization matrix C
   */
  base::OperationMatrix* computeRegularizationMatrix(Grid& grid);

  /**
   * Numerator of the system right hand side b contribution for 1st dataset
   */
  DataVector bNumP;

  /**
   * Denominator of the system right hand side b contribution for 1st dataset
   */
  DataVector bDenomP;

  /**
   * Numerator of the system right hand side b contribution for 2nd dataset
   */
  DataVector bNumQ;

  /**
   * Denominator of the system right hand side b contribution for 2nd dataset
   */
  DataVector bDenomQ;

  /**
   * Function that indicates whether a model is refinable at all (certain on/off settings do not
   * allow for refinement)
   * @return whether the model is refinable
   */
  bool isRefinable() override;
};
} /* namespace datadriven */
} /* namespace sgpp */
