/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ModelFittingDensityEstimation.hpp
 *
 * Created on: Jan 02, 2018
 *     Author: Kilian Röhner
 */

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>

using sgpp::base::DataMatrix;
using sgpp::base::Grid;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven {

// TODO(lettrich): allow different refinement techniques.
/**
 * Fitter object that encapsulates the usage of sparse grid density estimation with identity as
 * regularization.
 *
 * Allows usage of different grids, different solvers and different regularization techniques based
 * on the provided configuration objects.
 */
class ModelFittingDensityEstimation : public ModelFittingBaseSingleGrid {
 public:
  /**
   * Constructor
   *
   * @param config configuration object that specifies grid, refinement, and regularization
   */
  explicit ModelFittingDensityEstimation(const FitterConfigurationDensityEstimation& config);

  /**
   * Fit the grid to the given dataset by determining the weights of the initial grid by the
   * SGDE approach.
   * @param dataset the training dataset that is used to fit the model.
   */
  void fit(Dataset& dataset) override;

  /**
   * Improve accuracy of the fit on the given training data by adaptive refinement of the grid and
   * recalculate weights.
   * @return true if refinement could be performed based on the refinement configuration, else
   * false.
   */
  bool refine() override;

  void update(Dataset& dataset) override;

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

 private:
  /**
   * Count the amount of refinement operations performed on the current dataset.
   */
  size_t refinementsPerformed;

  /**
   * Reset the state of the object when a new dataset is used;
   */
  void resetState();

  // The offline object (contains decomposed matrix)
  std::unique_ptr<DBMatOffline> offline;

  // The online object
  std::unique_ptr<DBMatOnlineDE> online;
};
} /* namespace datadriven */
} /* namespace sgpp */
