// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/algorithm/DBMatObjectStore.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationClassification.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/globaldef.hpp>

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
   * @brief Constructor with specified object store.
   *
   * @param config Configuration object that specifies grid, refinement, and regularization
   * @param objectStore Offline object store for already decomposed offline objects.
   */
  explicit ModelFittingClassification(const FitterConfigurationClassification& config,
                                      std::shared_ptr<DBMatObjectStore> objectStore);

  /**
   * Fits the models for all classes based on the data given in the dataset parameter
   * @param dataset the training dataset that is used to fit the models
   */
  void fit(Dataset& dataset) override;
  void fit(Dataset& datasetP, Dataset& datasetQ) override {
    throw base::application_exception("This model requires a single input dataset");
  }

  /**
   * Improve the accuracy of the classification by refining or coarsening the grids of each class.
   * Coarsening is currently only implemented for RefinementFunctorType::Classification
   *  @return true if refinement could be performed for any grid based on the refinement
   * configuration, else false.
   */
  bool adapt() override;

  /**
   * Updates the models for each class based on new data (streaming or batch learning)
   * @param dataset the new data
   */
  void update(Dataset& dataset) override;
  void update(Dataset& datasetP, Dataset& datasetQ) override {
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
   * Should compute some kind of Residual to evaluate the fit of the model.
   *
   * In the case of density estimation, this is
   * || R * alpha_lambda - b_val ||_2
   *
   * This is useful for unsupervised learning models, where normal evaluation cannot be used as
   * there are no targets.
   *
   * For classification, this is not implemented, as accuracy should be used in this case.
   *
   * @param validationData Matrix for validation data
   *
   * @returns the residual score
   */
  double computeResidual(DataMatrix& validationData) const override {
    throw sgpp::base::not_implemented_exception(
        "ModelFittingDensityEstimationCombi::computeResidual() is not implemented!");
  }

  /**
   * Resets any trained representations of the model, but does not reset the entire state.
   *
   * Decompositions are not discarded, but can be reused.
   */
  void resetTraining() override;

  /**
   * Updates the regularization parameter lambda of the underlying model.
   *
   * @param lambda the new lambda parameter
   */
  void updateRegularization(double lambda) override;

  /**
   * Resets the state of the entire model
   */
  void reset() override;

  /**
   * store Fitter into text file in folder /datadriven/classificator/
   */
  void storeClassificator();

  /**
   * obtain the density estimation models per each class. To be used in VisualizerClassification
   */
  std::vector<std::unique_ptr<ModelFittingDensityEstimation>>* getModels();

  /**
   * obtain the index mapping for each label class. To be used in VisualizerClassification
   */
  std::map<double, size_t> getClassIdx();

#ifdef USE_SCALAPACK
  /**
   * @returns the BLACS process grid
   */
  std::shared_ptr<BlacsProcessGrid> getProcessGrid() const override;
#endif

 private:
  /**
   * @brief Offline object store for already decomposed offline objects.
   *
   */
  std::shared_ptr<DBMatObjectStore> objectStore;

  /**
   * @brief Flag to specify whether the instance has an object store.
   *
   */
  bool hasObjectStore;

  /**
   * Translates a class label to an index for the models vector. If the class is not present it will
   * create a new index for this class
   * @param label the label the translate
   * @return the index of this class label
   */
  size_t labelToIdx(double label);

  std::vector<double> getClassPriors() const;

  /**
   * Returns the refinement functor suitable for the model settings.
   * @param grids vector of pointers to grids for each class
   * @param surpluses vector of pointers to the surpluses for each class
   * @param priors vector of priors for each class
   * @return pointer to a refinement functor that suits the model settings
   */
  MultiGridRefinementFunctor* getRefinementFunctor(std::vector<Grid*> grids,
                                                   std::vector<DataVector*> surpluses,
                                                   std::vector<double> priors);

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
   * Initial size of the grids.
   */
  size_t initialGridSize;

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
