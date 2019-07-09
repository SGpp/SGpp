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

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDensityEstimation.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/ListNode.hpp>
#include <sgpp/base/tools/json/DictNode.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/algorithms/bhtsne/tsne.hpp>
#include <omp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>


using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using json::JSON;
using json::DictNode;
using json::ListNode;
using sgpp::datadriven::TSNE;


namespace sgpp {
namespace datadriven {


VisualizerDensityEstimation::VisualizerDensityEstimation(VisualizerConfiguration config) {
  this->config = config;
}

void VisualizerDensityEstimation::visualize(ModelFittingBase &model,
  size_t fold, size_t batch) {
  if (batch % config.getGeneralConfig().numBatches != 0) {
    return;
  }

  createOutputDirectory(fold, batch);

  omp_set_num_threads(static_cast<int> (config.getVisualizationParameters().numberCores));

  #pragma omp parallel sections
  {
    #pragma omp section
    {
    getLinearCuts(model);
    }

    #pragma omp section
    {
    getHeatmap(model);
    }

    #pragma omp section
    {
      if (config.getGeneralConfig().algorithm == "tsne") {
        runTsne(model);
      }
    }
    #pragma omp section
    {
      if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
        storeGrid(model);
      }
    }
  }
}

void VisualizerDensityEstimation::storeGrid(ModelFittingBase &model) {
  ModelFittingBaseSingleGrid* gridModel = dynamic_cast<ModelFittingBaseSingleGrid*>(&model);

  auto grid = gridModel->getGrid().clone();

  DataMatrix gridMatrix;

  grid->getStorage().getCoordinateArrays(gridMatrix);

  CSVTools::writeMatrixToCSVFile(config.getGeneralConfig().currentDirectory + "/grid", gridMatrix);
}

void VisualizerDensityEstimation::runTsne(ModelFittingBase &model) {
  DataMatrix data = model.getDataset()->getData();

  double* input = model.getDataset()->getData().data();

  DataVector evaluation(data.getNrows());

  model.evaluate(data, evaluation);

  size_t N = model.getDataset()->getNumberInstances();
  size_t D = model.getDataset()->getDimension();

  if ( D == 1 ) {
  std::cout << "The tsne algorithm can only be applied if "
  "the dimension is greater than 1" << std::endl;
  return;
  }
  double* output = reinterpret_cast<double*>(malloc(N * config.getVisualizationParameters().
    targetDimension * sizeof(double)));

  if ( D > config.getVisualizationParameters().targetDimension ) {
    std::cout << "Compressing with tsne to " <<
    std::to_string(config.getVisualizationParameters().targetDimension)
    << " dimensions" << std::endl;

    TSNE::run(input, N, D , output, config.getVisualizationParameters().targetDimension,
    config.getVisualizationParameters().perplexity, config.getVisualizationParameters().theta,
    config.getVisualizationParameters().seed, false,
    config.getVisualizationParameters().maxNumberIterations);
    D = config.getVisualizationParameters().targetDimension;
  } else {
    output = input;
  }

  DataMatrix compressedModel(output, N, D);

  compressedModel.appendCol(evaluation);

  if ( config.getGeneralConfig().targetFileType == VisualizationFileType::CSV ) {
    CSVTools::writeMatrixToCSVFile(config.getGeneralConfig().currentDirectory +
      "/tsneCompression", compressedModel);
  } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
    if (config.getVisualizationParameters().targetDimension != 2) {
    std::cout << "A json output is only available for compressions in 2 dimensions"
    "Storing the CSV instead" << std::endl;
    CSVTools::writeMatrixToCSVFile(config.getGeneralConfig().currentDirectory +
      "/tsneCompression", compressedModel);
    }
    storeTsneJson(compressedModel, model);
  }
}

void VisualizerDensityEstimation::getLinearCuts(ModelFittingBase &model) {
  std::cout << "Generating the linear cuts" << std::endl;

  auto nDimensions  = model.getFitterConfiguration().getGridConfig().dim_;
  DataMatrix cutMatrix(0, nDimensions);

  if ( nDimensions >= 2 ) {
    std::string command("mkdir ");
    command.append(config.getGeneralConfig().currentDirectory);
    command.append("/LinearCuts");
    system(command.data());

    if (nDimensions >=3) {
      getLinearCutsMore3D(cutMatrix, model);
    } else {
      getLinearCuts2D(cutMatrix, model);
    }
  } else {
    getLinearCuts1D(cutMatrix, model);
  }
}

void VisualizerDensityEstimation::getHeatmap(ModelFittingBase &model) {
  std::cout << "Generating the heatmaps" << std::endl;

  auto nDimensions = model.getFitterConfiguration().getGridConfig().dim_;

  if ( nDimensions == 1 ) {
    std::cout << "Heatmap generation is not available for models of 1 dimension" <<std::endl;
    return;
  }

  DataMatrix heatMapMatrix(0, nDimensions);
  if ( nDimensions >=3 ) {
    std::string command("mkdir ");
    command.append(config.getGeneralConfig().currentDirectory);
    command.append("/Heatmaps");
    system(command.data());

    if ( nDimensions >= 4 ) {
      getHeatmapMore4D(heatMapMatrix, model);
    } else if ( nDimensions == 3 ) {
      getHeatmap3D(heatMapMatrix, model);
    }
  } else {
    getHeatmap2D(heatMapMatrix, model);
  }
}

void VisualizerDensityEstimation::translateColumns(DataMatrix &matrix, size_t maxColumns) {
  DataMatrix temp(matrix);

  DataVector column(matrix.getNrows());

  for (size_t dimension=0; dimension < maxColumns-1; dimension++) {
    matrix.getColumn(dimension, column);

    temp.setColumn(dimension+1, column);
  }
  matrix.getColumn(maxColumns-1, column);
  temp.setColumn(0, column);
  matrix.copyFrom(temp);
}

void VisualizerDensityEstimation::translateColumnsRight(DataMatrix &matrix,
  std::vector <size_t> indexes) {
  DataMatrix temp(matrix);

  DataVector column(matrix.getNrows());

  for (size_t dimension = 0; dimension < indexes.size()-1; dimension++) {
  matrix.getColumn(indexes.at(dimension), column);
  temp.setColumn(indexes.at(dimension+1), column);
  }

  matrix.getColumn(indexes.back(), column);
  temp.setColumn(indexes.front(), column);
  matrix.copyFrom(temp);
}

void VisualizerDensityEstimation::translateColumnsLeft(DataMatrix &matrix,
  std::vector <size_t> indexes) {
  DataMatrix temp(matrix);

  DataVector column(matrix.getNrows());

  for (size_t dimension = indexes.size()-1; dimension > 0; dimension--) {
    DataVector column(matrix.getNrows());
    matrix.getColumn(indexes.at(dimension), column);
    temp.setColumn(indexes.at(dimension-1), column);
  }

  matrix.getColumn(indexes.front(), column);
  temp.setColumn(indexes.back(), column);
  matrix.copyFrom(temp);
}

void VisualizerDensityEstimation::updateIndexesCuts(std::vector <size_t> &indexes,
  DataMatrix &matrix) {
  if (indexes.at(2) < matrix.getNcols()-1) {
    indexes.at(2)++;
    swapColumns(matrix, indexes.at(2)-1, indexes.at(2));
  } else {
    if (indexes.at(1) < matrix.getNcols()-2) {
      indexes.at(1)++;
      indexes.at(2) = indexes.at(1)+1;

      swapColumns(matrix, indexes.at(2)-1, indexes.at(2));
      swapColumns(matrix, indexes.at(1)-1, indexes.at(1));

    } else {
      indexes.at(0)++;
      indexes.at(1) = indexes.at(0)+1;
      indexes.at(2) = indexes.at(1)+1;

      swapColumns(matrix, matrix.getNcols()-1, indexes.at(2));
      swapColumns(matrix, matrix.getNcols()-2, indexes.at(1));
      swapColumns(matrix, indexes.at(0)-1, indexes.at(0));
    }
  }
}

void VisualizerDensityEstimation::updateIndexesHeatmap(std::vector <size_t> &indexes,
  DataMatrix &matrix) {
  if (indexes.at(3) < matrix.getNcols()-1) {
  indexes.at(3)++;
  swapColumns(matrix, indexes.at(3)-1, indexes.at(3));
  } else {
    if (indexes.at(2) < matrix.getNcols()-2) {
    indexes.at(2)++;
    indexes.at(3) = indexes.at(2)+1;

    swapColumns(matrix, indexes.at(2)-1, indexes.at(2));
    swapColumns(matrix, matrix.getNcols()-1, indexes.at(3));
    } else {
        if (indexes.at(1) < matrix.getNcols()-3) {
          indexes.at(1)++;
          indexes.at(2) = indexes.at(1)+1;
          indexes.at(3) = indexes.at(2)+1;

          swapColumns(matrix, matrix.getNcols()-1, indexes.at(3));
          swapColumns(matrix, matrix.getNcols()-2, indexes.at(2));
          swapColumns(matrix, indexes.at(1)-1, indexes.at(1));
        } else {
          indexes.at(0)++;
          indexes.at(1) = indexes.at(0)+1;
          indexes.at(2) = indexes.at(1)+1;
          indexes.at(2) = indexes.at(1)+1;
          indexes.at(3) = indexes.at(2)+1;

          swapColumns(matrix, matrix.getNcols()-1, indexes.at(3));
          swapColumns(matrix, matrix.getNcols()-2, indexes.at(2));
          swapColumns(matrix, matrix.getNcols()-3, indexes.at(1));
          swapColumns(matrix, indexes.at(0)-1, indexes.at(0));
        }
    }
  }
}

void VisualizerDensityEstimation::swapColumns(DataMatrix &matrix, size_t col1, size_t col2) {
  DataVector temp1(matrix.getNrows());

  DataVector temp2(matrix.getNrows());

  matrix.getColumn(col1, temp1);
  matrix.getColumn(col2, temp2);
  matrix.setColumn(col2, temp1);
  matrix.setColumn(col1, temp2);
}

void VisualizerDensityEstimation::getLinearCutsMore3D(DataMatrix &cutMatrix,
  ModelFittingBase &model) {
  std::string outputDir(config.getGeneralConfig().currentDirectory + "/LinearCuts/");

  for (double dim1 = 0; dim1 <= 1; dim1 += 0.25) {
    for (double dim2 = 0; dim2 <= 1; dim2 += 0.25) {
      for (double dim3=0; dim3 <= 1; dim3 = round((dim3 + 0.01) * 100) / 100) {
        DataVector row(cutMatrix.getNcols(), 0.5);
        row.set(2, dim2);
        row.set(0, dim3);
        row.set(1, dim1);
        cutMatrix.appendRow(row);
      }
    }
  }

  std::vector <size_t> variableColumnIndexes = {0, 1, 2};

  while (variableColumnIndexes.at(2) < cutMatrix.getNcols()) {
    std::string subfolder(outputDir+"dimensions_"+
    std::to_string(variableColumnIndexes.at(0)+1)+"_"+
    std::to_string(variableColumnIndexes.at(1)+1)+"_"+
    std::to_string(variableColumnIndexes.at(2)+1));

    std::string command("mkdir " + subfolder);

    system(command.data());


    for (size_t combination = 0; combination < 3; combination++) {
      DataMatrix cutResults(cutMatrix);
      DataVector evaluation(cutMatrix.getNrows());

      model.evaluate(cutMatrix, evaluation);

      cutResults.appendCol(evaluation);

      translateColumnsRight(cutMatrix, variableColumnIndexes);
      if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
        CSVTools::writeMatrixToCSVFile(subfolder + "/Cut_var_dimension_" +
        std::to_string(variableColumnIndexes.at(combination)+1),
        cutResults);
      } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
        storeCutJson(cutResults, variableColumnIndexes, variableColumnIndexes.at(combination),
        subfolder + "/Cut_var_dimension_"+
        std::to_string(variableColumnIndexes.at(combination)+1));
      }
    }
  updateIndexesCuts(variableColumnIndexes, cutMatrix);
  }
}

void VisualizerDensityEstimation::getLinearCuts2D(DataMatrix &cutMatrix,
  ModelFittingBase &model) {
  std::string outputDir(config.getGeneralConfig().currentDirectory + "/LinearCuts/");

  for (double dim1 = 0; dim1 <= 1; dim1 += 0.25) {
    for (double dim2 = 0; dim2 <= 1; dim2 = round((dim2 + 0.01) * 100) / 100) {
        DataVector row(2, 0.5);
        row.set(0, dim2);
        row.set(1, dim1);
        cutMatrix.appendRow(row);
     }
  }

  std::vector <size_t> variableColumnIndexes = {0, 1, 2};

  for (size_t combination = 0; combination < 2; combination++) {
    DataMatrix cutResults(cutMatrix);
    DataVector evaluation(cutMatrix.getNrows());

    model.evaluate(cutMatrix, evaluation);

    cutResults.appendCol(evaluation);

    translateColumns(cutMatrix, cutMatrix.getNcols());

    if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
      CSVTools::writeMatrixToCSVFile(outputDir + "Cut_dimensions_1_2_variable_dimension" +
      std::to_string(combination+1),
      cutResults);
    } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
      storeCutJson(cutResults, variableColumnIndexes, variableColumnIndexes.at(combination),
      outputDir + "Cut_dimensions_1_2_variable_dimension" +
      std::to_string(combination+1));
    }
  }
}

void VisualizerDensityEstimation::getLinearCuts1D(DataMatrix &cutMatrix, ModelFittingBase &model) {
  std::string outputDir(config.getGeneralConfig().currentDirectory+"/");

  for (double dim1 = 0; dim1 <= 1; dim1 = round((dim1 + 0.01) * 100) / 100) {
    DataVector row(1, 0.5);
    row.set(0, dim1);
    cutMatrix.appendRow(row);
  }

  DataMatrix cutResults(cutMatrix);
  DataVector evaluation(cutMatrix.getNrows());

  model.evaluate(cutMatrix, evaluation);

  cutResults.appendCol(evaluation);
  if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
    CSVTools::writeMatrixToCSVFile(outputDir+"FittedModel", cutResults);
  } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
  storeCutJson(cutResults, outputDir+"FittedModel");
  }
}

void VisualizerDensityEstimation::getHeatmapMore4D(DataMatrix &heatMapMatrix,
ModelFittingBase &model) {
  std::string outputDir(config.getGeneralConfig().currentDirectory+"/Heatmaps/");
  for (double dim1 = 0; dim1 <= 1; dim1+=0.25) {
    for (double dim2 = 0; dim2 <= 1; dim2+=0.25) {
      for (double dim3 = 0; dim3 <= 1; dim3 = roundf((dim3+0.05) * 100) / 100) {
        for (double dim4 = 0; dim4 <= 1; dim4 = roundf((dim4+0.05) * 100) / 100) {
        DataVector row(heatMapMatrix.getNcols(), 0.5);
        row.set(0, dim3);
        row.set(1, dim4);
        row.set(2, dim1);
        row.set(3, dim2);

        heatMapMatrix.appendRow(row);
        }
      }
    }
  }

  std::vector <size_t> variableColumnIndexes = {0, 1, 2, 3};

  std::vector<size_t> workingIndexes = {0, 0, 0};

  while (variableColumnIndexes.at(3) < heatMapMatrix.getNcols()) {
    std::string subfolder(outputDir+"dimensions_"+
    std::to_string(variableColumnIndexes.at(0)+1)+"_"+
    std::to_string(variableColumnIndexes.at(1)+1)+"_"+
    std::to_string(variableColumnIndexes.at(2)+1)+"_"+
    std::to_string(variableColumnIndexes.at(3)+1));

    std::string command("mkdir "+subfolder);

    system(command.data());


    std::copy(variableColumnIndexes.begin()+1, variableColumnIndexes.end(),
    workingIndexes.begin());

    for (size_t iteration = 0; iteration < 2; iteration++) {
      for (size_t combination = 0; combination < 3; combination++) {
        DataMatrix heatMapResults(heatMapMatrix);
        DataVector evaluation(heatMapMatrix.getNrows());

        model.evaluate(heatMapMatrix, evaluation);

        heatMapResults.appendCol(evaluation);

        if (iteration == 0) {
          if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
            CSVTools::writeMatrixToCSVFile(subfolder+"/Heatmap_var_dimensions_"
            +std::to_string(variableColumnIndexes.at(0)+1)+"_"+
            std::to_string(variableColumnIndexes.at(combination+1)+1), heatMapResults);
          } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
            storeHeatmapJson(heatMapResults, model,
            variableColumnIndexes, variableColumnIndexes.at(0),
            variableColumnIndexes.at(combination+1),
            subfolder+"/Heatmap_var_dimensions_"
                       +std::to_string(variableColumnIndexes.at(0)+1)+"_"+
                       std::to_string(variableColumnIndexes.at(combination+1)+1));
          }
          translateColumnsRight(heatMapMatrix, workingIndexes);
        } else {
          if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
            CSVTools::writeMatrixToCSVFile(subfolder+"/Heatmap_var_dimensions_"
              +std::to_string(variableColumnIndexes.at((combination < 2)?1:2)+1)+"_"+
              std::to_string(variableColumnIndexes.at((combination < 1)?2:3)+1), heatMapResults);
          } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
            storeHeatmapJson(heatMapResults, model,
             variableColumnIndexes, variableColumnIndexes.at((combination < 2)?1:2),
             variableColumnIndexes.at((combination < 1)?2:3),
             subfolder+"/Heatmap_var_dimensions_"
             +std::to_string(variableColumnIndexes.at((combination < 2)?1:2)+1)+"_"+
             std::to_string(variableColumnIndexes.at((combination < 1)?2:3)+1));
          }
          translateColumnsLeft(heatMapMatrix, workingIndexes);
        }
      }

      if (iteration == 0) {
        translateColumnsRight(heatMapMatrix, variableColumnIndexes);
      }
    }
    translateColumnsLeft(heatMapMatrix, variableColumnIndexes);

    updateIndexesHeatmap(variableColumnIndexes, heatMapMatrix);
  }
}

void VisualizerDensityEstimation::getHeatmap3D(DataMatrix &heatMapMatrix, ModelFittingBase &model) {
  std::string outputDir(config.getGeneralConfig().currentDirectory+"/Heatmaps/");

  // Dummy to reutilize the storejson method
  std::vector <size_t> variableColumnIndexes = {0, 1, 2};

  for (double dim1 = 0; dim1 <= 1; dim1+=0.25) {
    for (double dim2 = 0; dim2 <= 1; dim2 = round((dim2+0.05) * 100) / 100) {
      for (double dim3 = 0; dim3 <= 1; dim3 = round((dim3+0.05) * 100) / 100) {
        DataVector row(3, 0.5);
        row.set(0, dim3);
        row.set(1, dim2);
        row.set(2, dim1);
        heatMapMatrix.appendRow(row);
      }
    }
  }

  for (size_t combination = 0; combination < 3; combination++) {
    DataMatrix heatMapResults(heatMapMatrix);
    DataVector evaluation(heatMapMatrix.getNrows());

    model.evaluate(heatMapMatrix, evaluation);

    heatMapResults.appendCol(evaluation);

    translateColumns(heatMapMatrix, heatMapMatrix.getNcols());
    if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
      CSVTools::writeMatrixToCSVFile(outputDir+"Heatmap_var_dimensions_"
      +std::to_string(combination+1)
      +"_"+((combination < 2)?std::to_string(combination+2):std::to_string(1)), heatMapResults);

    } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
      storeHeatmapJson(heatMapResults,
      model,
      variableColumnIndexes,
      variableColumnIndexes.at(combination),
      variableColumnIndexes.at((combination < 2)?combination+1:0),
      outputDir+"Heatmap_var_dimensions_"+std::to_string(combination+1)+"_"+
      ((combination < 2)?std::to_string(combination+2):std::to_string(1)));
    }
  }
}

void VisualizerDensityEstimation::getHeatmap2D(DataMatrix &heatMapMatrix,
  ModelFittingBase &model) {
  std::string outputDir(config.getGeneralConfig().currentDirectory+"/");
  for (double dim1 = 0; dim1 <= 1; dim1 = round((dim1+0.05) * 100) / 100) {
    for (double dim2 = 0; dim2 <= 1; dim2 = round((dim2+0.05) * 100) / 100) {
      DataVector row(2, 0.5);
      row.set(0, dim1);
      row.set(1, dim2);
      heatMapMatrix.appendRow(row);
    }
  }

  DataMatrix heatMapResults(heatMapMatrix);
  DataVector evaluation(heatMapMatrix.getNrows());

  model.evaluate(heatMapMatrix, evaluation);

  heatMapResults.appendCol(evaluation);

  if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
  CSVTools::writeMatrixToCSVFile(outputDir+"FittedModel", heatMapResults);
  } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
    storeHeatmapJson(heatMapResults,
    model, outputDir+"FittedModel");
  }
}



void VisualizerDensityEstimation::storeTsneJson(DataMatrix &matrix, ModelFittingBase &model) {
  JSON jsonOutput;

  jsonOutput.addListAttr("data");

  jsonOutput["data"].addDictValue();
  jsonOutput["data"][0].addIDAttr("type", "\"scatter\"");
  jsonOutput["data"][0].addIDAttr("mode", "\"markers\"");

  DataVector xCol(matrix.getNrows());

  matrix.getColumn(0, xCol);

  jsonOutput["data"][0].addIDAttr("x", xCol.toString());

  DataVector yCol(matrix.getNrows());

  matrix.getColumn(1, yCol);
  jsonOutput["data"][0].addIDAttr("y", yCol.toString());

  jsonOutput["data"][0].addDictAttr("marker");

  DataVector zCol(matrix.getNrows());

  matrix.getColumn(2, zCol);

  jsonOutput["data"][0]["marker"].addIDAttr("color", zCol.toString());

  jsonOutput["data"][0]["marker"].addIDAttr("colorscale", "\"Viridis\"");

  jsonOutput["data"][0]["marker"].addIDAttr("opacity", 0.8);

  jsonOutput["data"][0]["marker"].addIDAttr("showscale", true);

  jsonOutput["data"][0]["marker"].addDictAttr("colorbar");

  jsonOutput["data"][0]["marker"]["colorbar"].addDictAttr("title");
  jsonOutput["data"][0]["marker"]["colorbar"]["title"].addIDAttr("text", "\"Density value\"");

  jsonOutput.addDictAttr("layout");

  jsonOutput["layout"].addDictAttr("title");

  jsonOutput["layout"]["title"].addIDAttr("text", "\"TSNE Compression\"");

  jsonOutput.serialize(config.getGeneralConfig().currentDirectory + "/tsneCompression.json");
}

void VisualizerDensityEstimation::storeCutJson(DataMatrix &matrix, std::vector<size_t> indexes,
size_t &varDim, std::string filepath) {
  indexes.erase(std::find(indexes.begin(), indexes.end(), varDim));

  JSON jsonOutput;
  jsonOutput.addListAttr("data");

  jsonOutput.addDictAttr("layout");

  jsonOutput["layout"].addDictAttr("title");

  jsonOutput["layout"].addListAttr("annotations");

  jsonOutput["layout"]["title"].addIDAttr("text", "\"Linear Cuts: Variable dimension " +
  std::to_string(varDim + 1) + "\"");

  for (unsigned int graphNumber=0; graphNumber < matrix.getNrows()/101; graphNumber++) {
    DataMatrix temp(matrix.getNrows(), matrix.getNcols());

    temp.copyFrom(matrix);

    unsigned int beginIndex = graphNumber*101+1;

    temp.resizeToSubMatrix(beginIndex, 1, beginIndex + 100, matrix.getNcols());

    // Adding data trace
    jsonOutput["data"].addDictValue();

    jsonOutput["data"][graphNumber].addIDAttr("type", "\"scatter\"");
    jsonOutput["data"][graphNumber].addIDAttr("mode", "\"lines\"");

    std::string xAxis("x" + ((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    jsonOutput["data"][graphNumber].addIDAttr("xaxis",
      "\"x" + ((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");

    std::string yAxis("y"+((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    jsonOutput["data"][graphNumber].addIDAttr("yaxis",
      "\"y" + ((graphNumber == 0)?"":std::to_string(graphNumber + 1)) + "\"");

    DataVector xCol(temp.getNrows());

    temp.getColumn(varDim, xCol);

    jsonOutput["data"][graphNumber].addIDAttr("x", xCol.toString());

    DataVector yCol(temp.getNrows());

    temp.getColumn(temp.getNcols()-1, yCol);
    jsonOutput["data"][graphNumber].addIDAttr("y", yCol.toString());
    jsonOutput["data"][graphNumber].addIDAttr("showlegend", false);
    jsonOutput["data"][graphNumber].addIDAttr("hoverinfo", "\"x+y\"");


    // Layout part of the graph
    std::string xAxisName("xaxis" + ((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    std::string yAxisName("yaxis" + ((graphNumber == 0)?"":std::to_string(graphNumber+1)));

    jsonOutput["layout"].addDictAttr(xAxisName);
    jsonOutput["layout"][xAxisName].addIDAttr("anchor", "\""+yAxis+"\"");
    jsonOutput["layout"][xAxisName].addIDAttr("type", "\"linear\"");
    jsonOutput["layout"][xAxisName].addListAttr("domain");
    jsonOutput["layout"][xAxisName]["domain"].addIdValue(0.2*(graphNumber%5));
    jsonOutput["layout"][xAxisName]["domain"].addIdValue(0.2*((graphNumber%5)+1)-0.05);

    jsonOutput["layout"].addDictAttr(yAxisName);
    jsonOutput["layout"][yAxisName].addIDAttr("anchor", "\""+xAxis+"\"");
    jsonOutput["layout"][yAxisName].addIDAttr("type", "\"linear\"");

    if (matrix.getNcols() > 3) {
      jsonOutput["layout"][yAxisName].addListAttr("domain");
      jsonOutput["layout"][yAxisName]["domain"].addIdValue(0.2*(graphNumber/5));
      jsonOutput["layout"][yAxisName]["domain"].addIdValue(0.2*((graphNumber/5)+1)-0.1);
    }

    // Adding titles to subplots
    DataVector firstRow(temp.getNcols());
    temp.getRow(0, firstRow);

    // Adding titles to subplots
    std::string dim1Text(std::to_string(indexes.at(0) + 1));

    std::string dim1ValueText(std::to_string(firstRow.get(indexes.at(0))));
    dim1ValueText.erase(dim1ValueText.find_last_not_of('0') + 2, std::string::npos);

    dim1Text = "\"Dim " + dim1Text + "=" + dim1ValueText;

    std::string dim2Text = "";
    std::string dim2ValueText = "\"";

    if (matrix.getNcols() > 3) {
      dim2Text = std::to_string(indexes.at(1) + 1);
      dim2ValueText = std::to_string(firstRow.get(indexes.at(1)));
      dim2ValueText.erase(dim2ValueText.find_last_not_of('0') + 2, std::string::npos);
      dim2Text = ", Dim " + dim2Text + "= "+dim2ValueText + "\"";
    } else {
      dim2Text = dim2Text + dim2ValueText;
    }

    std::string subplot_title(dim1Text + dim2Text);

    jsonOutput["layout"]["annotations"].addDictValue();
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("x",
    (std::stod(jsonOutput["layout"][xAxisName]["domain"][0].get()) +
    std::stod(jsonOutput["layout"][xAxisName]["domain"][1].get()))/2);

    if (matrix.getNcols() > 3) {
      jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("y",
      (0.9-std::stod(jsonOutput["layout"][yAxisName]["domain"][0].get())));
    } else {
      jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("y", 1.0);
    }

    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("showarrow", false);
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("xanchor", "\"center\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("yanchor", "\"bottom\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("xref", "\"paper\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("yref", "\"paper\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("text", subplot_title);
  }
  std::cout << "Writing file " << filepath+".json" << std::endl;
  jsonOutput.serialize(filepath + ".json");
}

void VisualizerDensityEstimation::storeCutJson(DataMatrix &matrix, std::string filepath) {
  JSON jsonOutput;
  jsonOutput.addListAttr("data");

  jsonOutput.addDictAttr("layout");

  jsonOutput["layout"].addDictAttr("title");

  jsonOutput["layout"]["title"].addIDAttr("text", "\"Density Estimation: Fitted Function \"");

  // Adding data trace
  jsonOutput["data"].addDictValue();

  jsonOutput["data"][0].addIDAttr("type", "\"scatter\"");
  jsonOutput["data"][0].addIDAttr("mode", "\"lines\"");

  DataVector xCol(matrix.getNrows());

  matrix.getColumn(0, xCol);

  jsonOutput["data"][0].addIDAttr("x", xCol.toString());

  DataVector yCol(matrix.getNrows());

  matrix.getColumn(1, yCol);
  jsonOutput["data"][0].addIDAttr("y", yCol.toString());
  jsonOutput["data"][0].addIDAttr("showlegend", false);
  jsonOutput["data"][0].addIDAttr("hoverinfo", "\"x+y\"");

  std::cout << "Writing file " << filepath+".json" << std::endl;
  jsonOutput.serialize(filepath + ".json");
}

void VisualizerDensityEstimation::storeHeatmapJson(DataMatrix &matrix, ModelFittingBase &model,
std::vector<size_t> indexes, size_t &varDim1, size_t &varDim2, std::string filepath) {
  indexes.erase(std::find(indexes.begin(), indexes.end(), varDim1));

  indexes.erase(std::find(indexes.begin(), indexes.end(), varDim2));

  ModelFittingBaseSingleGrid* gridModel = dynamic_cast<ModelFittingBaseSingleGrid*>(&model);

  auto grid = gridModel->getGrid().clone();

  DataMatrix gridMatrix;

  grid->getStorage().getCoordinateArrays(gridMatrix);

  double maxValue = matrix.max(matrix.getNcols()-1);
  double minValue = matrix.min(matrix.getNcols()-1);

  JSON jsonOutput;
  jsonOutput.addListAttr("data");
  jsonOutput.addDictAttr("layout");
  jsonOutput["layout"].addDictAttr("title");

  if (gridMatrix.getNcols() >= 4) {
    jsonOutput["layout"].addIDAttr("height", "1500");
  }
  jsonOutput["layout"].addListAttr("annotations");
  jsonOutput["layout"]["title"].addIDAttr("text", "\"Heatmaps: Variable dimensions: " +
  std::to_string(varDim1 + 1) + " and " + std::to_string(varDim2 + 1) + "\"");

  for (unsigned int graphNumber = 0; graphNumber < matrix.getNrows()/441; graphNumber++) {
    DataMatrix temp(matrix.getNrows(), matrix.getNcols());

    temp.copyFrom(matrix);

    unsigned int beginIndex = graphNumber*441+1;

    temp.resizeToSubMatrix(beginIndex, 1, beginIndex+440, matrix.getNcols());

    // Adding data trace
    jsonOutput["data"].addDictValue();

    jsonOutput["data"][2*graphNumber].addIDAttr("type", "\"contour\"");

    std::string xAxis("x" + ((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    jsonOutput["data"][2*graphNumber].addIDAttr("xaxis",
    "\"x"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");

    std::string yAxis("y" + ((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    jsonOutput["data"][2*graphNumber].addIDAttr("yaxis",
    "\"y"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");

    DataVector xCol(temp.getNrows());

    temp.getColumn(varDim1, xCol);

    jsonOutput["data"][2*graphNumber].addIDAttr("x", xCol.toString());

    DataVector yCol(temp.getNrows());

    temp.getColumn(varDim2, yCol);

    jsonOutput["data"][2*graphNumber].addIDAttr("y", yCol.toString());

    DataVector zCol(temp.getNrows());

    temp.getColumn(temp.getNcols()-1, zCol);

    jsonOutput["data"][2*graphNumber].addIDAttr("z", zCol.toString());
    jsonOutput["data"][2*graphNumber].addIDAttr("showlegend", false);
    jsonOutput["data"][2*graphNumber].addIDAttr("hoverinfo", "\"x+y+z\"");

    jsonOutput["data"][2*graphNumber].addIDAttr("zmin", minValue);
    jsonOutput["data"][2*graphNumber].addIDAttr("zmax", maxValue);
    jsonOutput["data"][2*graphNumber].addIDAttr("colorscale", "\"Viridis\"");

    // Adding the grid
    DataVector firstRow(temp.getNcols());
    temp.getRow(0, firstRow);
    DataMatrix tempGrid(0, gridMatrix.getNcols());

    for (size_t index = 0; index < gridMatrix.getNrows(); index++) {
      DataVector row(gridMatrix.getNcols());
      gridMatrix.getRow(index, row);

      if (gridMatrix.getNcols() >= 4 ) {
        if (row.get(indexes.at(0)) == firstRow.get(indexes.at(0)) &
          row.get(indexes.at(1)) == firstRow.get(indexes.at(1))) {
          tempGrid.appendRow(row);
        }
      } else {
        if (row.get(indexes.at(0)) == firstRow.get(indexes.at(0))) {
          tempGrid.appendRow(row);
        }
      }
    }

    jsonOutput["data"].addDictValue();
    jsonOutput["data"][2*graphNumber+1].addIDAttr("type", "\"scatter\"");
    jsonOutput["data"][2*graphNumber+1].addIDAttr("mode", "\"markers\"");
    jsonOutput["data"][2*graphNumber+1].addIDAttr("xaxis",
      "\"x"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");
    jsonOutput["data"][2*graphNumber+1].addIDAttr("yaxis",
      "\"y"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");
    jsonOutput["data"][2*graphNumber+1].addDictAttr("marker");
    jsonOutput["data"][2*graphNumber+1]["marker"].addIDAttr("color", "\"red\"");

    DataVector xColGrid(tempGrid.getNrows());

    tempGrid.getColumn(varDim1, xColGrid);

    jsonOutput["data"][2*graphNumber+1].addIDAttr("x", xColGrid.toString());

    DataVector yColGrid(tempGrid.getNrows());

    tempGrid.getColumn(varDim2, yColGrid);

    jsonOutput["data"][2*graphNumber+1].addIDAttr("y", yColGrid.toString());
    jsonOutput["data"][2*graphNumber+1].addIDAttr("showlegend", false);
    jsonOutput["data"][2*graphNumber+1].addIDAttr("hoverinfo", "\"none\"");

    // Layout part of the graph
    std::string xAxisName("xaxis"+((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    std::string yAxisName("yaxis"+((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    jsonOutput["layout"].addDictAttr(xAxisName);
    jsonOutput["layout"][xAxisName].addIDAttr("anchor", "\"" + yAxis + "\"");
    jsonOutput["layout"][xAxisName].addIDAttr("type", "\"linear\"");
    jsonOutput["layout"][xAxisName].addListAttr("domain");
    jsonOutput["layout"][xAxisName]["domain"].addIdValue(0.2*(graphNumber%5));
    jsonOutput["layout"][xAxisName]["domain"].addIdValue(0.2*((graphNumber%5)+1)-0.05);


    jsonOutput["layout"].addDictAttr(yAxisName);
    jsonOutput["layout"][yAxisName].addIDAttr("anchor", "\""+xAxis+"\"");
    jsonOutput["layout"][yAxisName].addIDAttr("type", "\"linear\"");

    if (gridMatrix.getNcols() >= 4) {
      jsonOutput["layout"][yAxisName].addListAttr("domain");
      jsonOutput["layout"][yAxisName]["domain"].addIdValue(0.2*(graphNumber/5));
      jsonOutput["layout"][yAxisName]["domain"].addIdValue(0.2*((graphNumber/5)+1)-0.1);
    }

    // Adding titles to subplots
    std::string dim1Text(std::to_string(indexes.at(0)+1));

    std::string dim1ValueText(std::to_string(firstRow.get(indexes.at(0))));
    dim1ValueText.erase(dim1ValueText.find_last_not_of('0') + 2, std::string::npos);

    dim1Text = "\"Dim " +
    dim1Text + "=" + dim1ValueText;

    std::string dim2Text = "";
    std::string dim2ValueText = "\"";

    if (gridMatrix.getNcols() >= 4) {
      dim2Text = std::to_string(indexes.at(1) + 1);
      dim2ValueText = std::to_string(firstRow.get(indexes.at(1)));
      dim2ValueText.erase(dim2ValueText.find_last_not_of('0') + 2, std::string::npos);
      dim2Text = ", Dim " + dim2Text+"= " + dim2ValueText + "\"";
    } else {
      dim2Text = dim2Text + dim2ValueText;
    }

    std::string subplot_title(dim1Text + dim2Text);

    jsonOutput["layout"]["annotations"].addDictValue();
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("x",
    (std::stod(jsonOutput["layout"][xAxisName]["domain"][0].get()) +
    std::stod(jsonOutput["layout"][xAxisName]["domain"][1].get()))/2);

    if (gridMatrix.getNcols() >= 4) {
      jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("y",
      (0.9-std::stod(jsonOutput["layout"][yAxisName]["domain"][0].get())));
    } else {
      jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("y", 1.0);
    }

    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("showarrow", false);
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("xanchor", "\"center\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("yanchor", "\"bottom\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("xref", "\"paper\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("yref", "\"paper\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("text", subplot_title);
  }
  std::cout << "Writing file " << filepath + ".json" << std::endl;
  jsonOutput.serialize(filepath+".json");
}

void VisualizerDensityEstimation::storeHeatmapJson(DataMatrix &matrix, ModelFittingBase &model,
std::string filepath) {
  ModelFittingBaseSingleGrid* gridModel = dynamic_cast<ModelFittingBaseSingleGrid*>(&model);

  auto grid = gridModel->getGrid().clone();

  DataMatrix gridMatrix;

  grid->getStorage().getCoordinateArrays(gridMatrix);

  double maxValue = matrix.max();
  double minValue = matrix.min();

  JSON jsonOutput;
  jsonOutput.addListAttr("data");

  jsonOutput.addDictAttr("layout");

  jsonOutput["layout"].addDictAttr("title");

  jsonOutput["layout"].addListAttr("annotations");

  jsonOutput["layout"]["title"].addIDAttr("text", "\"Density Estimation: 2D Fitted Model\"");

  // Adding data trace
  jsonOutput["data"].addDictValue();

  jsonOutput["data"][0].addIDAttr("type", "\"contour\"");

  DataVector xCol(matrix.getNrows());

  matrix.getColumn(0, xCol);

  jsonOutput["data"][0].addIDAttr("x", xCol.toString());

  DataVector yCol(matrix.getNrows());

  matrix.getColumn(1, yCol);

  jsonOutput["data"][0].addIDAttr("y", yCol.toString());

  DataVector zCol(matrix.getNrows());

  matrix.getColumn(2, zCol);

  jsonOutput["data"][0].addIDAttr("z", zCol.toString());
  jsonOutput["data"][0].addIDAttr("showlegend", false);
  jsonOutput["data"][0].addIDAttr("hoverinfo", "\"x+y+z\"");
  jsonOutput["data"][0].addIDAttr("zmin", minValue);
  jsonOutput["data"][0].addIDAttr("zmax", maxValue);
  jsonOutput["data"][0].addIDAttr("colorscale", "\"Viridis\"");

  // Adding the grid
  jsonOutput["data"].addDictValue();

  jsonOutput["data"][1].addIDAttr("type", "\"scatter\"");

  jsonOutput["data"][1].addIDAttr("mode", "\"markers\"");

  jsonOutput["data"][1].addDictAttr("marker");

  jsonOutput["data"][1]["marker"].addIDAttr("color", "\"red\"");

  DataVector xColGrid(gridMatrix.getNrows());

  gridMatrix.getColumn(0, xColGrid);

  jsonOutput["data"][1].addIDAttr("x", xColGrid.toString());

  DataVector yColGrid(gridMatrix.getNrows());

  gridMatrix.getColumn(1, yColGrid);

  jsonOutput["data"][1].addIDAttr("y", yColGrid.toString());
  jsonOutput["data"][1].addIDAttr("showlegend", false);
  jsonOutput["data"][1].addIDAttr("hoverinfo", "\"none\"");

  std::cout << "Writing file " << filepath + ".json" << std::endl;
  jsonOutput.serialize(filepath + ".json");
}

}  // namespace datadriven
}  // namespace sgpp
