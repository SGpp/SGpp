// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_BOOST_GRAPH
#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationClustering.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/tools/Graph.hpp>

#include<list>
#include <queue>
#include <vector>

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
   * Runs the series of stpes to obtain the clustering model of the dataset based on the
   * algorithm using a density estimation model.
   * Requires only data samples and no targets (since those are irrelevant for the clustering)
   * @param newDataset the training dataset that is used to fit the model.
   */
  void fit(Dataset& newDataset) override;

  /**
   * Updates the model with new samples
   * @param newDataset Training datasetused to update the model
   */
  void update(Dataset& newDataset) override;

  double evaluate(const DataVector& sample) override;

  /**
   * Evaluates a series of samples and assign them a cluster according to the
   * classification model
   * @param samples Samples to evaluate
   * @param results Labels assigned to the samples
   */
  void evaluate(DataMatrix& samples, DataVector& results) override;

  /**
   * Performs a refinement given the new grid size and the points to coarsened
   * @return if the grid was refined (true)
   */
  bool refine() override;

  /**
   * Clears the model
   */
  void reset() override;

  /**
   * Resets any trained representations of the model, but does not reset the entire state.
   */
  void resetTraining() override;

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
  double computeResidual(DataMatrix &validationData) const override;

  /**
   * Updates the regularization parameter lambda of the underlying model.
   *
   * @param lambda the new lambda parameter
   */
  void updateRegularization(double lambda) override;

  /**
   * Returns the pointer to the pointer pointing to the density estimation model
   * @return Pointer of the unique pointer of the density estimation model
   */
  std::unique_ptr<ModelFittingDensityEstimation>* getDensityEstimationModel();

    /**
   * Returns the pointer to the pointer pointing to the classification estimation model
   * @return Pointer of the unique pointer of the classification estimation model
   */
  std::unique_ptr<ModelFittingClassification>* getClassificationModel();


  DataMatrix &getLabeledPoints();

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
   * VpTree to do nearest neighbors queries
   */
  std::unique_ptr<VpTree> vpTree;

  /**
   * Nearest Neighbors' Graph
   */
  std::shared_ptr<Graph> graph;

  /**
   * Graph that keeps a copy of the nearest neighbors graph and on which the pruning is done
   */
  std::shared_ptr<Graph> prunedGraph;

  /**
   * Classification Model
   */
  std::unique_ptr<ModelFittingClassification> classificationModel;

  /**
   * Method which updates the nearest neighbors graph.
   * @param newDataset New dataset to be added
   */
  void updateGraph(DataMatrix &newDataset);

  /**
   * Method that generates the density estimation model of the unlabeled dataset
   * @params dataset the training dataset that is used to fit the model.
   */
  void generateDensityEstimationModel(Dataset &dataset);

  /**
   * Method that generates the similarity graph of a given dataset based on the nearest
   * neighbors hyperparameter given in the configuration
   * @params dataset the training dataset that is used to generate the graph
   */
  void generateSimilarityGraph(Dataset &dataset);

  /**
   * Method that deletes from the graphs all of the points that do not have the minimun
   * density given by the user as an hyperparameter in the configuration
   * @params deletedNodes Matrix to store all of the indexes of the
   * deleted points for further labelling
   */
  void applyDensityThresholds(std::vector<size_t> &deletedNodes);

  /**
   * Method that detects all of the disconected components from the graph and assigns a label
   * to each of its corresponding points
   * @params deletedNodes List of nodes which were deleted from the graph
   */
  void detectComponentsAndLabel(std::vector<size_t> &deletedNodes);

  /**
   * Method that generates a classification model to cluster out of sample points based
   * on the cluster obtained by the training data.
   * @params labeledSamples The matrix containing the training data with labeles to generate
   * the classification model.
   */
  void generateClassificationModel(DataMatrix &labeledSamples);

  /**
   * Creates a density estimation model that fits the model settings.
   * @param densityEstimationConfig configuration for the density estimation
   * @return a new density estimation model
   */
  std::unique_ptr<ModelFittingDensityEstimation> createNewDensityModel(
      sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig);
};
}  // namespace datadriven
}  // namespace sgpp

#endif