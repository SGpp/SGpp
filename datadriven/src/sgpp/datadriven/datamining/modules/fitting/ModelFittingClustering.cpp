// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClustering.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCG.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCombi.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOffParallel.hpp>

#include <map>
#include <iostream>
#include <ctime>
#include <queue>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

ModelFittingClustering::ModelFittingClustering(
  const FitterConfigurationClustering& config)
  : ModelFittingBase(), refinementsPerformed{0} {
  this->verboseSolver = true;
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationClustering>(config));

  this-> densityEstimationModel = createNewDensityModel(
      dynamic_cast<sgpp::datadriven::FitterConfigurationDensityEstimation&>(*(this->config)));
  this->classificationModel = std::make_unique<ModelFittingClassification>(
      dynamic_cast<sgpp::datadriven::FitterConfigurationClassification&>(*(this->config)));
  this->graph = nullptr;
  this->vpTree = nullptr;
#ifdef USE_SCALAPACK
    auto& parallelConfig = this->config->getParallelConfig();
  if (parallelConfig.scalapackEnabled_) {
    processGrid = std::make_shared<BlacsProcessGrid>(config.getParallelConfig().processRows_,
                                                     config.getParallelConfig().processCols_);
  }
#endif
}

void ModelFittingClustering::fit(Dataset& newDataset) {
  reset();
  update(newDataset);
}

void ModelFittingClustering::update(Dataset& newDataset) {
  // Update the model
  dataset = &newDataset;
  generateDensityEstimationModel(newDataset);
  if (firstEpoch) {
    generateSimilarityGraph(newDataset);
  }
  std::vector<size_t> deletedNodes;
  applyDensityThresholds(deletedNodes);
  detectComponentsAndLabel(deletedNodes);
}

double ModelFittingClustering::evaluate(const DataVector &sample) {
  return 0.0;
}

void ModelFittingClustering::evaluate(DataMatrix& samples, DataVector& results) {
  return;
}

bool ModelFittingClustering::refine() {
  densityEstimationModel->refine();
  std::vector<size_t> deletedNodes;
  applyDensityThresholds(deletedNodes);
  detectComponentsAndLabel(deletedNodes);
  return true;
}

void ModelFittingClustering::reset() {
  densityEstimationModel->reset();
  classificationModel->reset();
  graph = nullptr;
}

std::unique_ptr<ModelFittingDensityEstimation> ModelFittingClustering::createNewDensityModel(
    sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig) {
  if (densityEstimationConfig.getGridConfig().generalType_ ==
      base::GeneralGridType::ComponentGrid) {
    return std::make_unique<ModelFittingDensityEstimationCombi>(densityEstimationConfig);
  }
  switch (densityEstimationConfig.getDensityEstimationConfig().type_) {
    case DensityEstimationType::CG: {
      return std::make_unique<ModelFittingDensityEstimationCG>(densityEstimationConfig);
    }
    case DensityEstimationType::Decomposition: {
#ifdef USE_SCALAPACK
      if (densityEstimationConfig.getParallelConfig().scalapackEnabled_) {
    return std::make_unique<ModelFittingDensityEstimationOnOffParallel>(densityEstimationConfig,
                                                                        processGrid);
  }
#endif  // USE_SCALAPACK
      return std::make_unique<ModelFittingDensityEstimationOnOff>(densityEstimationConfig);
    }
  }
  throw application_exception("Unknown density estimation type");
}

void ModelFittingClustering::generateDensityEstimationModel(Dataset &dataset) {
  densityEstimationModel->update(dataset);
}

void ModelFittingClustering::generateSimilarityGraph(Dataset &dataset) {
  if (graph == nullptr) {
    this->vpTree = std::make_unique<VpTree>(dataset.getData());
    graph = std::make_shared<Graph>(dataset.getData().getNrows());
    DataVector currrentRow(dataset.getData().getNcols());

    for (size_t index = 0; index < vpTree->getStoredItems().getNrows(); index++) {
      vpTree->getStoredItems().getRow(index, currrentRow);
      auto nearestNeighbors = vpTree->getNearestNeighbors(currrentRow,
          config->getClusteringConfig().noNearestNeighbors);
      graph->createEdges(index, nearestNeighbors);
    }
    std::cout << "Num of vertices: " << boost::num_vertices(*(graph->getGraph())) << std::endl;
    std::cout << "Num of edges " << boost::num_edges(*(graph->getGraph())) << std::endl;

  } else {
    updateGraph(dataset.getData());
  }
}

void ModelFittingClustering::updateGraph(DataMatrix &newDataset) {
  clock_t start = std::clock();
  vpTree->update(newDataset);
  DataVector currrentRow(newDataset.getNcols());
  graph = std::make_shared<Graph>(vpTree->getStoredItems().getNrows());

  for (size_t index = 0; index < vpTree->getStoredItems().getNrows(); index++) {
    vpTree->getStoredItems().getRow(index, currrentRow);
    auto nearestNeighbors = vpTree->getNearestNeighbors(currrentRow,
        config->getClusteringConfig().noNearestNeighbors);
    graph->createEdges(index, nearestNeighbors);
  }
  clock_t end = std::clock();
  std::cout << "Vertices in *graph: " << boost::num_vertices(*graph->getGraph()) << std::endl;
  std::cout << "Edges in *graph: " << boost::num_edges(*graph->getGraph()) << std::endl;
  std::cout << "Graph updated the 'hard' way in  " <<
            std::to_string(static_cast<double>(end - start) / CLOCKS_PER_SEC) << " seconds"
            << std::endl;
}


void ModelFittingClustering::applyDensityThresholds(std::vector<size_t> &deletedNodes) {
  prunedGraph = std::make_shared<Graph>(*graph);
  DataVector evaluation(vpTree->getStoredItems().getNrows());

  densityEstimationModel->evaluate((vpTree->getStoredItems()), evaluation);
  evaluation.normalize();

  for (size_t index = 0; index < evaluation.size(); index++) {
    if (evaluation.get(index) < config->getClusteringConfig().densityThreshold) {
      prunedGraph->removeVertex(index-deletedNodes.size());
      deletedNodes.push_back(index);
    }
  }

  std::cout << "Number of deleted vertices: " << deletedNodes.size() <<std::endl;
  std::cout << "Remaining vertices: " <<
  boost::num_vertices(*(prunedGraph->getGraph())) <<std::endl;
}

void ModelFittingClustering::detectComponentsAndLabel(std::vector<size_t> &deletedNodes) {
  std::map<UndirectedGraph::vertex_descriptor, size_t> componentsMap;
  auto numberComponents = prunedGraph->getConnectedComponents(componentsMap);
  std::cout << "Number of found components: "<< numberComponents <<std::endl;

  std::map<int, size_t> numbersPerCluster;

  Dataset* classificationDataSet = new Dataset(vpTree->getStoredItems().getNrows(),
      vpTree->getStoredItems().getNcols());

  DataMatrix& samples = classificationDataSet->getData();
  DataVector& labels = classificationDataSet->getTargets();

  samples.copyFrom(vpTree->getStoredItems());

  for (size_t index = 0, graphIndex = 0; index < vpTree->getStoredItems().getNrows(); index++) {
    if (std::find(deletedNodes.begin(), deletedNodes.end(), index) != deletedNodes.end()) {
      labels.set(index, static_cast<double>(-1));
      numbersPerCluster[-1]++;
    } else {
      auto vertex_descriptor = boost::vertex(graphIndex++, *(prunedGraph->getGraph()));
      labels.set(index, static_cast<double>(componentsMap[vertex_descriptor]));
      numbersPerCluster[static_cast<int>(componentsMap[vertex_descriptor])]++;
    }
  }

  for (auto&cluster : numbersPerCluster) {
    std::cout << "Label " << cluster.first << ": " << cluster.second <<std::endl;
  }

  classificationModel->fit(*classificationDataSet);
}

std::unique_ptr<ModelFittingDensityEstimation>*
    ModelFittingClustering::getDensityEstimationModel() {
  return &densityEstimationModel;
}

std::unique_ptr<ModelFittingClassification>*
    ModelFittingClustering::getClassificationModel() {
  return &classificationModel;
}

DataMatrix &ModelFittingClustering::getLabeledPoints() {
  return vpTree->getStoredItems();
}
}  // namespace datadriven
}  // namespace sgpp
