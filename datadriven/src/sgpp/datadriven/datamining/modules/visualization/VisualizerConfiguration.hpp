// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/visualization/VisualizationParameters.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizationGeneralConfig.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

class VisualizerConfiguration{
 public:
  /**
   * Default Constructor
   */
  VisualizerConfiguration() = default;

  /**
   * Default Constructor
   */
  ~VisualizerConfiguration() = default;

  /**
   * set default values for all members based on the desired scenario.
   */
  void setupDefaults();

  /**
   * obtain parameters from a parser
   * @param parser: the parser object to read from
   */
  void readParams(const DataMiningConfigParser &parser);

  /**
   * read general configuration parameters
   */
  VisualizationGeneralConfig &getGeneralConfig();

  /**
   * read general configuration parameters
   */
  VisualizationParameters &getVisualizationParameters();

 protected:
  /**
   * Contains general configuration for the Visualizer
   */
  VisualizationGeneralConfig generalConfig;

  /**
   *  Contains the numerical parameters used to run the visualization
   *  algorithm
   */
  VisualizationParameters visualizationParameters;
};

}  // namespace datadriven
}  // namespace sgpp
