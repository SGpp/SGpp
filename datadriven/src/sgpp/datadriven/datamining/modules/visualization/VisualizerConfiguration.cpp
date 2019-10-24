// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizationGeneralConfig.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizationParameters.hpp>

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>


namespace sgpp {
namespace datadriven {

void VisualizerConfiguration::setupDefaults() {
  srand(static_cast<unsigned int>(time(NULL)));

  std::cout << "Setting up defaults parameters for Visualizer" << std::endl;
  generalConfig.algorithm = std::vector<std::string>({"tsne", "heatmaps", "linearcuts"});
  generalConfig.targetDirectory = "./output";
  generalConfig.targetFileType = VisualizationFileType::CSV;
  generalConfig.numBatches = 1;

  visualizationParameters.perplexity = 30;
  visualizationParameters.theta = 0.5;
  visualizationParameters.maxNumberIterations = 1000;
  visualizationParameters.seed = rand();
  visualizationParameters.targetDimension = 2;
  visualizationParameters.numberCores = 1;
  std::cout << "Setting up defaults done" << std::endl;
}

void VisualizerConfiguration::readParams(const DataMiningConfigParser &parser) {
  setupDefaults();
  parser.getVisualizationGeneralConfig(generalConfig, generalConfig);
  parser.getVisualizationParameters(visualizationParameters, visualizationParameters);
}

VisualizationGeneralConfig &VisualizerConfiguration::getGeneralConfig() {
  return generalConfig;
}

VisualizationParameters &VisualizerConfiguration::getVisualizationParameters() {
  return visualizationParameters;
}

}  // namespace datadriven
}  // namespace sgpp
