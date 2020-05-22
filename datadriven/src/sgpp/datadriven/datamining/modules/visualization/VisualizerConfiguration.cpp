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
  srand(static_cast<unsigned int>(time(nullptr)));

  std::cout << "Setting up defaults parameters for Visualizer" << std::endl;
  generalConfig.algorithm_ = std::vector<std::string>({"tsne", "heatmaps", "linearcuts"});
  generalConfig.targetDirectory_ = "./output";
  generalConfig.targetFileType_ = VisualizationFileType::CSV;
  generalConfig.numBatches_ = 1;

  visualizationParameters.perplexity_ = 30;
  visualizationParameters.theta_ = 0.5;
  visualizationParameters.maxNumberIterations_ = 1000;
  visualizationParameters.seed_ = rand();
  visualizationParameters.targetDimension_ = 2;
  visualizationParameters.numberCores_ = 1;
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
