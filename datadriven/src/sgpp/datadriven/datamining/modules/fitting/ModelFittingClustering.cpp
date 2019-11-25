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
  : refinementsPerformed{0} {
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
  //labeledSamples.resize(vpTree->getStoredItems().getNrows(), vpTree->getStoredItems().getNcols());
  //detectComponentsAndLabel(deletedNodes);
  //std::cout<<labeledSamples.toString();
}

double ModelFittingClustering::evaluate(const DataVector &sample) {
  return 0.0;
}

void ModelFittingClustering::evaluate(DataMatrix& samples, DataVector& results) {
  return;
}

bool ModelFittingClustering::refine() {
 // densityEstimationModel->refine();
  //std::vector<size_t> deletedNodes;
  //applyDensityThresholds(deletedNodes);
  //labeledSamples.resize(vpTree->getStoredItems().getNrows(), vpTree->getStoredItems().getNcols());
  //detectComponentsAndLabel(deletedNodes);

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

    for (size_t index = 0; index < dataset.getData().getNrows(); index++) {
      dataset.getData().getRow(index, currrentRow);
      auto nearestNeighbors = vpTree->getNearestNeighbors(currrentRow,
          config->getClusteringConfig().noNearestNeighbors);
      currentNearestNeighbors.push_back(nearestNeighbors);
      graph->createEdges(index, nearestNeighbors);
    }
    std::cout << "Num of vertices: " << boost::num_vertices(*(graph->getGraph()))<<std::endl;
    std::cout << "Num of edges " << boost::num_edges(*(graph->getGraph()))<<std::endl;

  } else {
    updateGraphOnline(dataset.getData());
    //updateHard(dataset.getData());
  }
}

void ModelFittingClustering::updateHard(DataMatrix &newDataset) {
  clock_t start = std::clock();
  vpTree->updateHard(newDataset);
  currentNearestNeighbors.clear();
  DataVector currrentRow(newDataset.getNcols());
  graph = std::make_shared<Graph>(vpTree->getStoredItems().getNrows());

  for (size_t index = 0; index<vpTree->getStoredItems().getNrows(); index++) {
    vpTree->getStoredItems().getRow(index, currrentRow);
    auto nearestNeighbors = vpTree->getNearestNeighbors(currrentRow,
        config->getClusteringConfig().noNearestNeighbors);
    currentNearestNeighbors.push_back(nearestNeighbors);
    graph->createEdges(index, nearestNeighbors);
  }
  clock_t end = std::clock();
  std::cout << "Vertices in *graph: " << boost::num_vertices(*graph->getGraph())<<std::endl;
  std::cout << "Edges in *graph: " << boost::num_edges(*graph->getGraph())<<std::endl;
  std::cout << "Graph updated the 'hard' way in  " <<
            std::to_string(static_cast<double>(end - start) / CLOCKS_PER_SEC) << " seconds"
            << std::endl;
}

void ModelFittingClustering::updateGraphOnline(DataMatrix &newDataset) {
  clock_t start = std::clock();
  DataVector currrentRow(newDataset.getNcols());
  DataVector currentHeapVector(newDataset.getNcols());

  size_t lastIndex = currentNearestNeighbors.size();
  vpTree->update(newDataset);
  std::cout<<"Adding vertexes online"<<std::endl;
  for (size_t index = 0; index<newDataset.getNrows(); index++) {
    graph->addVertex();
  }

  for (size_t index = 0; index<newDataset.getNrows(); index++) {
    newDataset.getRow(index, currrentRow);
    auto nearestNeighbors = vpTree->getNearestNeighbors(currrentRow,
        config->getClusteringConfig().noNearestNeighbors);
    currentNearestNeighbors.push_back(nearestNeighbors);
    graph->createEdges(lastIndex+index, nearestNeighbors);
    for (size_t indexHeaps = 0; indexHeaps < lastIndex; indexHeaps++ ) {
      vpTree->getStoredItems().getRow(currentNearestNeighbors[indexHeaps].top().index, currentHeapVector);
      double distance = vpTree->euclideanDistance(currrentRow, currentHeapVector);
      if (distance < currentNearestNeighbors[indexHeaps].top().distance) {
        graph->deleteEdge(indexHeaps, currentNearestNeighbors[indexHeaps].top().index);
        currentNearestNeighbors[indexHeaps].pop();
        currentNearestNeighbors[indexHeaps].push(
            VpHeapItem(lastIndex+index, distance));
        graph->addEdge(indexHeaps, lastIndex+index);
      }
    }
  }

  std::cout << "Vertices in graph: " << boost::num_vertices(*(graph->getGraph()))<<std::endl;
  std::cout << "Edges in graph: " << boost::num_edges(*(graph->getGraph()))<<std::endl;
  clock_t end = std::clock();
  std::cout << "Graph updated in " <<
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
      prunedGraph->disconnectVertex(index-deletedNodes.size());
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

  vpTree->getStoredItems() = classificationDataSet->getData();
  DataVector& labels = classificationDataSet->getTargets();

  for (size_t index = 0, graphIndex = 0; index < vpTree->getStoredItems().getNrows(); index++) {
    //size_t index = boost::get(boost::vertex_index, *(prunedGraph->getGraph()), p.first);
    if(std::find(deletedNodes.begin(), deletedNodes.end(), index) != deletedNodes.end()) {
      // labels for deleted nodes
      labels.set(index, -1);
      numbersPerCluster[-1]++;
    } else {
      auto vertex_descriptor = boost::vertex(graphIndex++, *(prunedGraph->getGraph()));
      labels.set(index, componentsMap[vertex_descriptor]);
      numbersPerCluster[componentsMap[vertex_descriptor]]++;
    }
  }

  for(auto&cluster: numbersPerCluster) {
    std::cout << "Label "<<cluster.first << ": "<< cluster.second <<std::endl;
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

}  // namespace datadriven
}  // namespace sgpp
