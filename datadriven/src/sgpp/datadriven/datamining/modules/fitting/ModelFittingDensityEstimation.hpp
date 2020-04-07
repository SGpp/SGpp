// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>
#include <sgpp/globaldef.hpp>

#include <list>
#include <memory>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Abstract super class to encapsulate density estimation models such as using offline/-online
 * splitting or conjugate gradients in order to solve the system.
 */
class ModelFittingDensityEstimation : public ModelFittingBaseSingleGrid {
 public:
  /**
   * Default constructor
   */
  ModelFittingDensityEstimation();

  /**
   * Fit the grid to the given dataset by determining the surpluses of the initial grid by the SGDE
   * approach. Requires only data samples and no targets (since those are irrelevant for the density
   * estimation whatsoever)
   * @param dataset the training dataset that is used to fit the model.
   */
  virtual void fit(DataMatrix& dataset) = 0;
  virtual void fit(DataMatrix& datasetP, DataMatrix& datasetQ) = 0;

  void fit(Dataset& dataset) override = 0;
  void fit(Dataset& datasetP, Dataset& datasetQ) override = 0;

  /**
   * Updates the model based on new data samples (streaming, batch learning). Requires only the data
   * samples and no targets (since those are irrelevant for the density estimation whatsoever)
   * @param samples the new data samples
   */
  virtual void update(DataMatrix& samples) = 0;
  virtual void update(DataMatrix& samplesP, DataMatrix& samplesQ) = 0;

  void update(Dataset& dataset) override = 0;
  void update(Dataset& datasetP, Dataset& datasetQ) override = 0;

  double evaluate(const DataVector& sample) override = 0;

  void evaluate(DataMatrix& samples, DataVector& results) override = 0;

  /**
   * Performs refinement and coarsening given the new grid size and the points to coarsened
   * @param newNoPoints the grid size after refinement and coarsening
   * @param deletedGridPoints a list of indexes for grid points that will be removed
   * @return if the grid was refined (true)
   */
  virtual bool adapt(size_t newNoPoints, std::vector<size_t>& deletedGridPoints) = 0;

  /**
   * Improve accuracy of the fit on the given training data by adaptive refinement or coarsening of
   * the grid and recalculate weights.
   * @return true if refinement or coarsening could be performed based on the refinement
   * configuration, else false.
   */
  bool adapt() override;

  /**
   * Computes a residual to evaluate the fit of the model.
   *
   * This is useful for density estimation, because other Scores cannot be used as there are no
   * targets.
   *
   * @param validationData Matrix for validation data
   */
  double computeResidual(DataMatrix& validationData) const override = 0;

  /**
   * Updates the regularization parameter lambda of the underlying model.
   *
   * @param lambda the new lambda parameter
   */
  void updateRegularization(double lambda) override = 0;

  /**
   * Resets any trained representations of the model, but does not reset the entire state.
   *
   * The decomposition from the offline phase is kept so that it can be reused.
   */
  void resetTraining() override = 0;

 protected:
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
   * Function that indicates whether a model is refinable at all (certain on/off settings do not
   * allow for refinement)
   * @return whether the model is refinable
   */
  virtual bool isRefinable() = 0;

  /**
   * Count the amount of refinement operations performed on the current dataset.
   */
  size_t refinementsPerformed;

  /**
   * Initial number of grid points
   */
  size_t initialGridSize;
};
} /* namespace datadriven */
} /* namespace sgpp */
