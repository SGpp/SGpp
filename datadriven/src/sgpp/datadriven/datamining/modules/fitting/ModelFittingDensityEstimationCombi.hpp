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
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Abstract super class to encapsulate density estimation models working with multiple Grids.
 */

class ModelFittingDensityEstimationCombi : public ModelFittingDensityEstimation {
 public:
  /**
   * Default constructor
   */
  ModelFittingDensityEstimationCombi();

  /**
   * Constructor from a FitterConfig
   * @param config FitterConfig
   */
  explicit ModelFittingDensityEstimationCombi(FitterConfigurationDensityEstimation& config);

  /**
   * Fit the grids to the given dataset by determining the weights of the initial grid by the
   * SGDE approach.
   * @param newDataset the training dataset that is used to fit the model.
   */
  void fit(Dataset& newDataset);

  /**
   * Fit the grids to the given dataset by determining the surpluses of the initial grid by the
   * SGDE approach. Requires only data samples and no targets (since those are irrelevant for the
   * density estimation whatsoever)
   * @param newDataset the training dataset that is used to fit the model.
   */
  void fit(DataMatrix& newDataset);

  void fit(Dataset&, Dataset&) {
    throw base::application_exception("This model requires a single input dataset");
  }
  void fit(DataMatrix&, DataMatrix&) {
    throw base::application_exception("This model requires a single input dataset");
  }

  void update(Dataset& dataset);

  /**
   * Updates the model based on new data samples (streaming, batch learning). Requires only
   * the data samples and no targets (since those are irrelevant for the density estimation
   * whatsoever)
   * @param samples the new data samples
   */
  void update(DataMatrix& samples);

  void update(Dataset&, Dataset&) {
    throw base::application_exception("This model requires a single input dataset");
  }
  void update(DataMatrix&, DataMatrix&) {
    throw base::application_exception("This model requires a single input dataset");
  }

  /**
   * Evaluate the fitted density at a single data point - requires a trained grid.
   * @param sample vector with the coordinates in all dimensions of that sample.
   * @return evaluation of the trained grid.
   */
  double evaluate(const DataVector& sample);

  /**
   * Evaluate the fitted density on a set of data points - requires a trained grid.
   * @param samples matrix where each row represents a sample and the columns contain the
   * coordinates in all dimensions of that sample.
   * @param results vector where each row will contain the evaluation of the respective sample on
   * the current model.
   */
  void evaluate(DataMatrix& samples, DataVector& results);

  /**
   * Refines the component with the biggest error
   * @return if a component was refined
   */
  bool refine();

  /**
   * Useless at the moment
   */
  bool refine(size_t newNoPoints, std::list<size_t>* deletedGridPoints);

  /**
   * Resets the state of the entire model
   */
  void reset() override;

 protected:
  /**
   * Contains the component grids witch form the sparse grids
   */
  std::vector<std::unique_ptr<ModelFittingDensityEstimation>> components;

  /**
   * Contains all level vector and weights of the current component grid set
   */
  std::vector<combiConfig> componentConfigs;

  /**
   * Contains the status of the component grids.
   * true: fitted
   * false: unfitted
   */
  std::vector<bool> fitted;

  /**
   * Delivers the initial level vectors and weighs and manages refinements
   */
  CombiConfigurator configurator;

  bool isRefinable();

  /**
   * Creates a density estimation model that fits the model settings.
   * @param densityEstimationConfig configuration for the density estimation
   * @return a new density estimation model
   */
  std::unique_ptr<ModelFittingDensityEstimation> createNewModel(
      sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig);

  void addNewModel(const combiConfig combiconfig);

  /**
   * @param indexRev the reverse index (distance from the end of the vector) of the component that
   * must be removed
   */
  void removeModel(size_t indexRev);

  // sgpp::DataMatrix datamatrix;
};

}  // namespace datadriven
}  // namespace sgpp
