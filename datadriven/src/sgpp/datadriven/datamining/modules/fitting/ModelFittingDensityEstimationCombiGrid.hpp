/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ModelFittingDensityEstimationMultiGrid.hpp
 *
 *  Created on: Jan 15, 2019
 *      Author: nico
 */

#pragma once

#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <sgpp/globaldef.hpp>

#include <list>

namespace sgpp {
namespace datadriven {

/**
 * Abstract super class to encapsulate density estimation models working with multiple Grids.
 */

class ModelFittingDensityEstimationCombiGrid : public ModelFittingDensityEstimation {
 public:
  /**
   * Default constructor
   */
  ModelFittingDensityEstimationCombiGrid();

  /**
   * Constructor from a FitterConfig
   * @param config FitterConfig
   */
  ModelFittingDensityEstimationCombiGrid(FitterConfigurationDensityEstimation& config);

  /*
   * need to implement it here because DesityEstimation works with DataMatrix and ModelFittingBase
   * with Dataset,
   * may override virtual fit function in modelfittingbasesinglegrid
   */
  /**
   * Fit the grids to the given dataset by determining the surpluses of the initial grid by the
   * SGDE approach. Requires only data samples and no targets (since those are irrelevant for the
   * density estimation whatsoever)
   * @param dataset the training dataset that is used to fit the model.
   */
  void fit(DataMatrix& dataset);

  void fit(Dataset& dataset);

  bool refine();

  void reset() override;

  bool refine(size_t newNoPoints, std::list<size_t>* deletedGridPoints);

  double evaluate(const DataVector& sample);

  void evaluate(DataMatrix& samples, DataVector& results);

  /**
   * Updates the model based on new data samples (streaming, batch learning). Requires only
   * the data samples and no targets (since those are irrelevant for the density estimation
   * whatsoever)
   * @param samples the new data samples
   */
  void update(DataMatrix& samples);

  void update(Dataset& dataset);

 protected:
  std::vector<std::unique_ptr<ModelFittingDensityEstimationOnOff>> models;

  std::vector<double> weights;

  bool isRefinable();
};

}  // namespace datadriven
}  // namespace sgpp
