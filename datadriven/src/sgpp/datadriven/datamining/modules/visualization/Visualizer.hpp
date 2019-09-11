/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Visualizer.hpp
 *
 *  Created on: 16th Jun 2019
 *      Author: Vincent Bautista
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>

#include "string"

namespace sgpp {
namespace datadriven {

class Visualizer{
 public:
  Visualizer();
 /**
  * Virtual destructor
  */
  virtual ~Visualizer() = default;

  /**
   * Method to run the visualization process for a given batch and fold
   * @param model The model used to evaluate the visualization
   * @param dataSource The datasource from where the data points are obtained
   * @param fold The current fold being processed
   * @param batch The current batch being processed
   */
  virtual void runVisualization(ModelFittingBase &model, DataSource &dataSource,
    size_t fold, size_t batch) = 0;

  /**
   * Get the configuration of the visualizer object.
   * @return configuration of the visualizer object
   */
  const VisualizerConfiguration &getVisualizerConfiguration() const;


 protected:
  /**
   * Method to create a folder based on the operating system of the user
   * @param folder_path absolute or relative path of the folder to be created
   */
  void createFolder(std::string folder_path);
  /**
   * Configuration object for the fitter.
   */
  VisualizerConfiguration config;

  /**
   * Matrix to store the data from the datasource
   */
  DataMatrix originalData;
  /**
   * Matrix to store the data provided by the tsne algorithm
   */
  DataMatrix tsneCompressedData;

  /**
   * Resolution used in the graphs
   */
  size_t resolution;

  /**
   * Variable to store the current folder in which the visualizer is working on
   */
  std::string currentDirectory;
};

}  // namespace datadriven
}  // namespace sgpp
