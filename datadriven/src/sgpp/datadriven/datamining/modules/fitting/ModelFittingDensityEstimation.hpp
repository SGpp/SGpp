// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>

#include <list>

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
   * Fit the grid to the given dataset by determining the surpluses of the initial grid by the
   * SGDE approach. Requires only data samples and no targets (since those are irrelevant for the
   * density estimation whatsoever)
   * @param dataset the training dataset that is used to fit the model.
   */
  virtual void fit(DataMatrix& dataset) = 0;

  void fit(Dataset& dataset) override = 0;

  /**
   * Updates the model based on new data samples (streaming, batch learning). Requires only
   * the data samples and no targets (since those are irrelevant for the density estimation
   * whatsoever)
   * @param samples the new data samples
   */
  virtual void update(DataMatrix& samples) = 0;

  void update(Dataset& dataset) override = 0;

  double evaluate(const DataVector& sample) override = 0;

  void evaluate(DataMatrix& samples, DataVector& results) override = 0;

  /**
   * Performs a refinement given the new grid size and the points to coarsened
   * @param newNoPoints the grid size after refinement and coarsening
   * @param deletedGridPoints a list of indexes for grid points that will be removed
   * @return if the grid was refined (true)
   */
  virtual bool refine(size_t newNoPoints, std::list<size_t>* deletedGridPoints) = 0;

  /**
   * Improve accuracy of the fit on the given training data by adaptive refinement of the grid and
   * recalculate weights.
   * @return true if refinement could be performed based on the refinement configuration, else
   * false.
   */
  bool refine() override;

  /**
   * Computes a residual to evaluate the fit of the model.
   *
   * This is useful for density estimation, because other Scores cannot be used as
   * there are no targets.
   *
   * @param validationData Matrix for validation data
   */
  double computeResidual(DataMatrix& validationData) const override = 0;

  /**
   * Updates the regularization parameter lambda of the underlying model.
   *
   * @param double the new lambda parameter
   */
  void updateRegularization(double lambda) override = 0;

  /**
   * Returns the refinement functor suitable for the model settings.
   * @return pointer to a refinement functor that suits the model settings
   */
  sgpp::base::RefinementFunctor* getRefinementFunctor();

 protected:
  /**
   * Count the amount of refinement operations performed on the current dataset.
   */
  size_t refinementsPerformed;

  /**
   * Function that indicates whether a model is refinable at all (certain on/off settings do not
   * allow for refinement)
   * @return whether the model is refinable
   */
  virtual bool isRefinable() = 0;
};
} /* namespace datadriven */
} /* namespace sgpp */
