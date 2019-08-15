/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * VisualizerDensityEstimation.hpp
 *
 *  Created on: 16th Jun 2019
 *      Author: Vincent Bautista
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/Visualizer.hpp>
#include <sgpp/datadriven/tools/CSVTools.hpp>
#include <vector>
#include <string>

namespace sgpp {
namespace datadriven {

class VisualizerDensityEstimation:public Visualizer {
 public:

  VisualizerDensityEstimation()= default;
  /**
   * Constructor given a configuration
   * @param config. The VisualizerConfiguration object which contains
   * the configuration to run the visualization module
   */
  explicit VisualizerDensityEstimation(VisualizerConfiguration config);

  ~VisualizerDensityEstimation() = default;

  void runVisualization(ModelFittingBase &model, DataSource &dataSource,
    size_t fold, size_t batch) override;

 protected:
  void runTsne(ModelFittingBase &model, DataSource &dataSource,
    size_t fold, size_t batch) override;

  void getHeatmap(ModelFittingBase &model, std::string currentDirectory);

  void getLinearCuts(ModelFittingBase &model, std::string currentDirectory);

  void storeGrid(ModelFittingBase &model, std::string currentDirectory);

  void getLinearCutsMore3D(ModelFittingBase &model, std::string currentDirectory);
  void getLinearCuts1D(ModelFittingBase &model, std::string currentDirectory);
  void getLinearCuts2D(ModelFittingBase &model, std::string currentDirectory);

  void getHeatmapMore4D(ModelFittingBase &model, std::string currentDirectory);
  void getHeatmap3D(ModelFittingBase &model, std::string currentDirectory);
  void getHeatmap2D(ModelFittingBase &model, std::string currentDirectory);


  void translateColumns(DataMatrix &matrix, size_t maxColumns);
  void translateColumnsRight(DataMatrix &matrix, std::vector<size_t> indexes);
  void translateColumnsLeft(DataMatrix &matrix, std::vector<size_t> indexes);
  void updateIndexesCuts(std::vector<size_t> &columnIndexes, DataMatrix &matrix);
  void updateIndexesHeatmap(std::vector<size_t> &columnIndexes, DataMatrix &matrix);
  void swapColumns(DataMatrix &matrix, size_t col1, size_t col2);

  void storeTsneJson(DataMatrix &matrix, ModelFittingBase &model, std::string currentDirectory);

  void storeCutJson(DataMatrix &matrix,
    std::vector<size_t> indexes, size_t &varDim, std::string filepath);
  void storeCutJson(DataMatrix &matrix, std::string filepath);
  void storeHeatmapJson(DataMatrix &matrix, ModelFittingBase &model,
    std::vector<size_t> indexes, size_t &varDim1, size_t &varDim2, std::string filepath);
  void storeHeatmapJson(DataMatrix &matrix, ModelFittingBase &model, std::string filepath);

  /**
   * Method which builds the matrices used to generate the cuts and the
   * heatmaps
   */
  void initializeMatrices(ModelFittingBase &model);

   /**
    * Variable to store the cut matrix to be evaluated
    */
   DataMatrix cutMatrix;

   /**
    * Variable to store the heatmap matrix to be evaluated
    */
   DataMatrix heatMapMatrix;
};

}  // namespace datadriven
}  // namespace sgpp
