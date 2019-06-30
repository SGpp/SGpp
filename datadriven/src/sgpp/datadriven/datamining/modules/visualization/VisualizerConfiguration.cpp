/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * VisualizerConfiguration.hpp
 *
 *  Created on: 16th Jun 2019
 *      Author: Vincent Bautista
 */



#include <string>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>

#include <sgpp/datadriven/datamining/modules/visualization/VisualizationGeneralConfig.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizationParameters.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {


 void VisualizerConfiguration::setupDefaults(){

  std::cout <<"Setting up defaults parameters for Visualizer"<<std::endl;
  generalConfig.algorithm= "tsne";
  generalConfig.targetFile="./output";
  generalConfig.targetFileType = VisualizationFileType::CSV;

  visualizationParameters.perplexity=30;
  visualizationParameters.theta=0.5;
  visualizationParameters.maxNunmberIterations=1000;
  visualizationParameters.seed=30;
  visualizationParameters.targetDimension=2;
  visualizationParameters.numberCores=1;
  std::cout <<"Setting up defaults done"<<std::endl;
 }

 void VisualizerConfiguration::readParams(const DataMiningConfigParser &parser) {
   setupDefaults();

   parser.getVisualizationGeneralConfig(generalConfig, generalConfig);

   parser.getVisualizationParameters(visualizationParameters,visualizationParameters);
 }

 VisualizationGeneralConfig &VisualizerConfiguration::getGeneralConfig()
 {
   return generalConfig;
 }

 VisualizationParameters &VisualizerConfiguration::getVisualizationParameters()
 {
   return visualizationParameters;
 }

}
}
