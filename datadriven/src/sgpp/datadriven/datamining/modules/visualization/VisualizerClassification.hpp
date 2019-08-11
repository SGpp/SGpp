/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * VisualizerDensityEstimation.hpp
 *
 *  Created on: 18th July 2019
 *      Author: Vincent Bautista
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDensityEstimation.hpp>
#include <sgpp/datadriven/tools/CSVTools.hpp>
#include <vector>
#include <string>

namespace sgpp {
namespace datadriven {

class VisualizerClassification:public VisualizerDensityEstimation {
 public:
  VisualizerClassification()= default;
   /**
    * Constructor given a configuration
    * @param config. The VisualizerConfiguration object which contains
    * the configuration to run the visualization module
    */
   explicit VisualizerClassification(VisualizerConfiguration config);

   ~VisualizerClassification() = default;

   void runVisualization(ModelFittingBase &model, DataSource &dataSource,
     size_t fold, size_t batch) override;

 protected:
  void getHeatmapsClassification(ModelFittingBase &model);

  void getHeatmapMore4DClassification(DataMatrix &heatMapMatrix,
   ModelFittingBase &model);

  void getHeatmap2DClassification(DataMatrix &heatMapMatrix,
     ModelFittingBase &model);

  void storeTsneJson(DataMatrix &matrix, ModelFittingBase &model);

  void storeHeatmapJsonClassification(DataMatrix &matrix, ModelFittingBase &model,
  std::string filepath);

  void storeHeatmapJsonClassification(DataMatrix &matrix, ModelFittingBase &model,
  std::vector<size_t> indexes, size_t &varDim1, size_t &varDim2, std::string filepath);

 private:
  /**
   * Method which returns a vector with the values of all the classes used in the model
   */
  void getClasses(DataVector &column, size_t numClasses, DataVector &classes);

  /**
   * List of colors to add to the grids per class in high dimensional cases
   */
  std::vector<std::string> colors = {"red", "darkviolet", "orange", "palegreen",
                    "plum", "purple", "chocolate", "darkcyan", "gold","tomato"};
};

}  // namespace datadriven
}  // namespace sgpp
