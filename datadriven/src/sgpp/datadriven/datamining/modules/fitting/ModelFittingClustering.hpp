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
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/tools/Graph.hpp>

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
   * @param newDataset the training dataset that is used to fit the model.
   */
  void fit(Dataset& newDataset) override;

  void update(Dataset& newDataset) override;

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
   * VpTree to do neareast neigbors queries
   */
  std::unique_ptr<VpTree> vpTree;

  std::vector<std::priority_queue<VpHeapItem>> currentNearestNeighbors;
  /**
   * Nearest Neighbors' Graph
   */

  std::shared_ptr<Graph> graph;

    /**
   * Graph structure to store the nearest neighbors graph
   */
   std::shared_ptr<Graph> prunedGraph;

  /**
   * Classification Model
   */
  std::unique_ptr<ModelFittingClassification> classificationModel;

  /**
   * Matrix to store the current labeled samples
   */
  sgpp::base::DataMatrix labeledSamples;

  void updateGraphOnline(DataMatrix &newDataset);

  void updateHard(DataMatrix &newDataset);

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
   * density given by the user as an hyperparameter inthe configuration
   * @params deletedPoints Matrix to store all of the deleted points for further labelling
   */
  void applyDensityThresholds(std::vector<size_t> &deletedNodes);

  /**
   * Method that detects all of the disconected components from the graph and assigns a label
   * to each of its corresponding points
   */
  void detectComponentsAndLabel(std::vector<size_t> &deletedNodes);

  /**
   * Method that generates a classification model to cluster out of sample points based
   * on the cluster obtained by the training data.
   * @params labeledSamples The matrix containing the training data with labeles to generate
   * the classification model.
   */
  void generateClassificationModel(DataMatrix &labeledSamples);

  std::unique_ptr<ModelFittingDensityEstimation> createNewDensityModel(
      sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig);

};
}  // namespace datadriven
}  // namespace sgpp
