// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationClustering.hpp>

#include<list>

namespace sgpp {
namespace datadriven {

/**
 * Fitter object that encapsulates density based classification using instances of
 * ModelFittingDensityEstimation for each class.
 */
class ModelFittingClustering : public ModelFittingBase {
 public:
  /**
   * Constructor
   *
   * @param config configuration object that specifies grid, refinement, and regularization
   */
  explicit ModelFittingClustering(const FitterConfigurationClustering& config);

  /**
   * Runs the series of stepes to obtain the clustering of the dataset based on the
   * algorithm using a density estimation model.
   * Requires only data samples and no targets (since those are irrelevant for the clustering)
   * Labeles for each cluster are generated during the fitting
   * @param dataset the training dataset that is used to fit the model.
   */
  void fit(Dataset& dataset) override;

  void update(Dataset& dataset) override;

  double evaluate(const DataVector& sample) override;

  void evaluate(DataMatrix& samples, DataVector& results) override;

  /**
   * Performs a refinement given the new grid size and the points to coarsened
   * @return if the grid was refined (true)
   */
  bool refine() override ;

  void reset() override;

  std::unique_ptr<ModelFittingDensityEstimation>* getDensityEstimationModel();

  std::unique_ptr<ModelFittingClassification>* getClassificationModel();

 protected:
  /**
   * Count the amount of refinement operations performed on the current dataset.
   */
  size_t refinementsPerformed;

 private:
  /**
  * Density Estimation Model
  */
  std::unique_ptr<ModelFittingDensityEstimation> densityEstimationModel;

  /**
   * Classification Model
   */
  std::unique_ptr<ModelFittingClassification> classificationModel;

  /**
   * Method that generates the density estimation model of the unlabeled dataset
   * @params dataset the training dataset that is used to fit the model.
   */
  void generateDensityEstimationModel(Dataset dataset);

  /**
   * Method that generates the similarity graph of a given dataset based on the nearest
   * neighbors hyperparameter given in the configuration
   * @params dataset the training dataset that is used to generate the graph
   */
  void generateSimilarityGraph(Dataset dataset);

  /**
   * Method that deletes from the graphs all of the points that do not have the minimun
   * density given by the user as an hyperparameter inthe configuration
   * @params deletedPoints Matrix to store all of the deleted points for further labelling
   */
  void applyDensityThresholds(DataMatrix &deletedPoints);

  /**
   * Method that detects all of the disconected components from the graph and assigns a label
   * to each of its corresponding points
   * @params labeledSamples DataMatrix containing the training dataset and
   * the labels of the clustering
   */
  void detectComponents(DataMatrix &labeledSamples);

  /**
   * Method that generates a classification model to cluster out of sample points based
   * on the cluster obtained by the training data.
   * @params labeledSamples The matrix containing the training data with labeles to generate
   * the classification model.
   */
  void generateClassificationModel(DataMatrix &labeledSamples);
};
}  // namespace datadriven
}  // namespace sgpp
