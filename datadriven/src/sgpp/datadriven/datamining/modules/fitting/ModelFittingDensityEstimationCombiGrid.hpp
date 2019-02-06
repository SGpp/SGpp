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
#include <sgpp/datadriven/datamining/configuration/CombiConfigurator.hpp>
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

  /**
   * Fit the grids to the given dataset by determining the weights of the initial grid by the
   * SGDE approach.
   * @param dataset the training dataset that is used to fit the model.
   */
  void fit(Dataset& newDataset);

  /**
   * Fit the grids to the given dataset by determining the surpluses of the initial grid by the
   * SGDE approach. Requires only data samples and no targets (since those are irrelevant for the
   * density estimation whatsoever)
   * @param dataset the training dataset that is used to fit the model.
   */
  void fit(DataMatrix& NewDataset);

  bool refine();

  bool refine(size_t newNoPoints, std::list<size_t>* deletedGridPoints);

  void update(Dataset& dataset);

  /**
   * Updates the model based on new data samples (streaming, batch learning). Requires only
   * the data samples and no targets (since those are irrelevant for the density estimation
   * whatsoever)
   * @param samples the new data samples
   */
  void update(DataMatrix& samples);

  double evaluate(const DataVector& sample);

  void evaluate(DataMatrix& samples, DataVector& results);

  void reset() override;

  /**
   * Creates a density estimation model that fits the model settings.
   * @param densityEstimationConfig configuration for the density estimation
   * @return a new density estimation model
   */
  std::unique_ptr<ModelFittingDensityEstimation> createNewModel(
      sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig);

  void addNewModel(combiConfig config);

 protected:
  std::vector<std::unique_ptr<ModelFittingDensityEstimation>> models;

  std::vector<double> weights;

  bool isRefinable();
};

}  // namespace datadriven
}  // namespace sgpp
