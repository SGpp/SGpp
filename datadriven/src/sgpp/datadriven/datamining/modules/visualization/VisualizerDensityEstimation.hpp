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
   * @param config The VisualizerConfiguration object which contains
   * the configuration to run the visualization module
   */
  explicit VisualizerDensityEstimation(VisualizerConfiguration config);

  /**
   * Default destructor
   */
  ~VisualizerDensityEstimation() = default;

  /**
   * Method to run the visualization process for a given batch and fold
   * @param model The model used to evaluate the visualization
   * @param dataSource The datasource from where the data points are obtained
   * @param fold The current fold being processed
   * @param batch The current batch being processed
   */
  void runVisualization(ModelFittingBase &model, DataSource &dataSource,
    size_t fold, size_t batch) override;

 protected:
  /*
   * Method which starts the heatmap generation for Density Estimation Models
   * @model The model used to evaluate the heatmap
   * @currentDirectory The current directory to store the heatmap results
   */
  void getHeatmap(ModelFittingBase &model, std::string currentDirectory);

  /*
   * Method which starts the linear cut generation for Density Estimation Models
   * @model The model used to evaluate the linear cuts
   * @currentDirectory The current directory to store the linear cuts results
   */
  void getLinearCuts(ModelFittingBase &model, std::string currentDirectory);

  /*
   * Method which stores the coordinates of the grid points
   * of a Density Estimation model
   * @model The model from where the grid is obtained
   * @currentDirectory The current directory to store the grid data
   */
  void storeGrid(ModelFittingBase &model, std::string currentDirectory);

  /*
   * Method which generates the linear cuts graphs for models of 3 or more dimensions
   * @model the model used to evaluate the linear cuts
   * @currentDirectory The current directory to store the linear cuts results
   */
  void getLinearCutsMore3D(ModelFittingBase &model, std::string currentDirectory);

  /*
   * Method which generates the linear cuts graphs for models of 1 dimension
   * @model the model used to evaluate the linear cuts
   * @currentDirectory The current directory to store the linear cuts results
   */
  void getLinearCuts1D(ModelFittingBase &model, std::string currentDirectory);

  /*
   * Method which generates the linear cuts graphs for models of 2 dimensions
   * @model the model used to evaluate the linear cuts
   * @currentDirectory The current directory to store the linear cuts results
   */
  void getLinearCuts2D(ModelFittingBase &model, std::string currentDirectory);

  /*
   * Method which generates the heatmap of models of 4 or more dimensions
   * @model The model used to evaluate the heatmap
   * @currentDirectory The current directory to store the heatmap results
   */
  void getHeatmapMore4D(ModelFittingBase &model, std::string currentDirectory);

  /*
   * Method which generates the heatmap of models of 3 dimensions
   * @model The model used to evaluate the heatmap
   * @currentDirectory The current directory to store the heatmap results
   */
  void getHeatmap3D(ModelFittingBase &model, std::string currentDirectory);

  /*
   * Method which generates the heatmap of models of 2 dimensions
   * @model The model used to evaluate the heatmap
   * @currentDirectory The current directory to store the heatmap results
   */
  void getHeatmap2D(ModelFittingBase &model, std::string currentDirectory);

  /*
   * Method which shifts one position the columns of a matrix from left to right
   * in a circular fashion until the column given by the parameters maxColumns
   * @param matrix The matrix to be shifted
   * @param maxColumn The max number of columns used when shifting
   */
  void translateColumns(DataMatrix &matrix, size_t maxColumns);

  /*
   * Method which shifts the columns given by the vector indexes
   * of a matrix from left to right
   * in a circular fashion. If indexes are <1,3,6> Then column 1 will be shifted to
   * 3, 3 to 6 and 6 to 1.
   * @param matrix The matrix to be shifted
   * @param indexes Vector containing the columns to shift
   */
  void translateColumnsRight(DataMatrix &matrix, std::vector<size_t> indexes);

  /*
   * Method which shifts the columns given by the vector indexes
   * of a matrix from right to left
   * in a circular fashion. If indexes are <1,3,6> Then column 1 will be shifted to
   * 6, 6 to 3 and 1 to 6.
   * @param matrix The matrix to be shifted
   * @param indexes Vector containing the columns to shift
   */
  void translateColumnsLeft(DataMatrix &matrix, std::vector<size_t> indexes);

  /*
   * Method to update the columns indexes to be shifted when generating
   * the linear cuts
   * @param columnIndexes Vector to update
   * @param matrix The matrix whose columns are being shifted
   */
  void updateIndexesCuts(std::vector<size_t> &columnIndexes, DataMatrix &matrix);

  /*
   * Method to update the columns indexes to be shifted when generating
   * the heatmap
   * @param columnIndexes Vector to update
   * @param matrix The matrix whose columns are being shifted
   */
  void updateIndexesHeatmap(std::vector<size_t> &columnIndexes, DataMatrix &matrix);

  /**
   * Method to swap to columns of a matrix
   * @param matrix The matrix whose columns are to be swaped
   * @param col1 The column identifier to be swaped with the column identified by col2
   * @param col2 The colum identifier to be swaped with the column identified by col1
   */
  void swapColumns(DataMatrix &matrix, size_t col1, size_t col2);

  /**
   * Method to generate and store in json  format for the
   * plotly library the output of the tsne algorithm
   * @param matrix Matrix with the content to be stored
   * @param model Model used in the evaluation
   * @param currentDirectory The current directory to store the json file
   */
  void storeTsneJson(DataMatrix &matrix, ModelFittingBase &model, std::string currentDirectory);

  /**
   * Method to generate and store in json  format for the
   * plotly library the output of the linear cuts for models of 2 or more dimensions
   * @param matrix Matrix with the content to be stored
   * @param indexes Vectors containing the dimensions used when generating these cuts
   * @param varDim THe dimension number varying and whose evaluation is shown in the model
   * @param filepath The current directory to store the json file
   */
  void storeCutJson(DataMatrix &matrix,
    std::vector<size_t> indexes, size_t &varDim, std::string filepath);

  /**
   * Method to generate and store in json  format for the
   * plotly library the output of the linear cuts for models of 1 dimension
   * @param matrix Matrix with the content to be stored
   * @param filepath The current directory to store the json file
   */
  void storeCutJson(DataMatrix &matrix, std::string filepath);

  /**
   * Method to generate and store in json  format for the
   * plotly library the output of the hetamaps for models of 3 or more dimensions
   * @param matrix Matrix with the content to be stored
   * @param model The model used when evaluating the heatmaps
   * @param indexes Vectors containing the dimensions used when generating these heatmaps
   * @param varDim1 The first dimension number varying and whose evaluation
   * is shown in the model
   * @param varDim2 The second dimension number varying and whose evaluation
   * is shown in the model
   * @param filepath The current directory to store the json file
   */
  void storeHeatmapJson(DataMatrix &matrix, ModelFittingBase &model,
    std::vector<size_t> indexes, size_t &varDim1, size_t &varDim2, std::string filepath);

  /**
   * Method to generate and store in json  format for the
   * plotly library the output of the heatmaps for models of 2 dimensions
   * @param matrix Matrix with the content to be stored
   * @param model The model used when evaluating the heatmaps
   * @param filepath The current directory to store the json file
   */
  void storeHeatmapJson(DataMatrix &matrix, ModelFittingBase &model, std::string filepath);

  /**
   * Method which builds the matrices used to generate the cuts and the
   * heatmaps
   * @param model The model used to evaluate the linear cuts and the heatmaps
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
