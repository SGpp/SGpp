/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ModelFittingClassification.hpp
 *
 *  Created on: Jul 1, 2018
 *      Author: dominik
 */

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationClassification.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>


#include <map>
#include <memory>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;

namespace sgpp {
namespace datadriven {

/**
 * Fitter object that encapsulates density based classification using instances of
 * ModelFittingDensityEstimation for each class.
 */
class ModelFittingClassification : public ModelFittingBase {
 public:
  /**
   * Constructor
   *
   * @param config configuration object that specifies grid, refinement, and regularization
   */
  explicit ModelFittingClassification(const FitterConfigurationClassification& config);

  /**
   * Fits the models for all classes based on the data given in the dataset parameter
   * @param dataset the training dataset that is used to fit the models
   */
  void fit(Dataset& dataset) override;
  void fit(Dataset&, Dataset&) override {
    throw base::application_exception("This model requires a single input dataset");
  }

  /**
   * Improve the accuracy of the classification by refining the grids of each class
   * @return true if refinement could be performed for any grid based on the refinement
   * configuration, else false.
   */
  bool refine() override;

  /**
   * Updates the models for each class based on new data (streaming or batch learning)
   * @param dataset the new data
   */
  void update(Dataset& dataset) override;
  void update(Dataset&, Dataset&) override {
    throw base::application_exception("This model requires a single input dataset");
  }

  /**
   * Predict the class of a data sample based on the density of the sample for each model
   * @param sample the sample point to classify
   * @return the predicted class label
   */
  double evaluate(const DataVector& sample) override;

  /**
   * Predicts the class for a set of data points based on the learned densities for each class
   * @param samples matrix where each row represents a data sample
   * @param results vector to output the predicted classes
   */
  void evaluate(DataMatrix& samples, DataVector& results) override;

  /**
   * Resets the state of the entire model
   */
  void reset() override;

  /*
   * store Fitter into text file in folder /datadriven/classificator/
   */
  void storeClassificator();

#ifdef USE_SCALAPACK
    /**
   * @returns the BLACS process grid
   */
  std::shared_ptr<BlacsProcessGrid> getProcessGrid() const override;
#endif

 private:
  /**
   * Translates a class label to an index for the models vector. If the class is not present
   * it will create a new index for this class
   * @param label the label the translate
   * @return the index of this class label
   */
  size_t labelToIdx(double label);

  /**
   * Returns the refinement functor suitable for the model settings.
   * @param grids vector of pointers to grids for each class
   * @param surpluses vector of pointers to the suprluses for each class
   * @return pointer to a refinement functor that suits the model settings
   */
  MultiGridRefinementFunctor* getRefinementFunctor(std::vector<Grid*> grids,
                                                   std::vector<DataVector*> surpluses);

  /**
   * Creates a density estimation model that fits the model settings.
   * @param densityEstimationConfig configuration for the density estimation
   * @return a new density estimation model
   */
  std::unique_ptr<ModelFittingDensityEstimation> createNewModel(
      sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig);

  /**
   * Count the amount of refinement operations performed on the current dataset.
   */
  size_t refinementsPerformed;

  /**
   * Models for each class
   */
  std::vector<std::unique_ptr<ModelFittingDensityEstimation>> models;

  /**
   * Index map for each class
   */
  std::map<double, size_t> classIdx;

  /**
   * Number of instances for each class (used for prior)
   */
  std::vector<size_t> classNumberInstances;

#ifdef USE_SCALAPACK
  /**
   * BLACS process grid for ScaLAPACK version
   */
  std::shared_ptr<BlacsProcessGrid> processGrid;
#endif
};
} /* namespace datadriven */
} /* namespace sgpp */
