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

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerClassification.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/ListNode.hpp>
#include <sgpp/base/tools/json/DictNode.hpp>
#include <omp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <math.h>

namespace sgpp {
namespace datadriven {
VisualizerClassification::VisualizerClassification(VisualizerConfiguration config) {
  this->config = config;
}


void VisualizerClassification::runVisualization(ModelFittingBase &model, DataSource &dataSource,
  size_t fold, size_t batch) {
  if (batch % config.getGeneralConfig().numBatches != 0) {
    return;
  }

  if(fold==0 && batch == 0) {
   originalData = dataSource.getAllSamples()->getData();
   resolution = pow(2,model.getFitterConfiguration().getGridConfig().level_+2);
  }

  ModelFittingClassification* classificationModel = dynamic_cast<ModelFittingClassification*>(&model);
  auto models = classificationModel->getModels();

  createOutputDirectory(fold, batch);

  omp_set_num_threads(static_cast<int> (config.getVisualizationParameters().numberCores));

  #pragma omp parallel sections
  {
    #pragma omp section
    {
     getHeatmapsClassification(model);
    }
    #pragma omp section
    {
       if (config.getGeneralConfig().algorithm == "tsne") {
         if(fold == 0 && batch == 0) {
           runTsne(model, dataSource, fold, batch);
         }
         if (originalData.getNcols() >= 1 ) {
           DataVector evaluation(originalData.getNrows());
           model.evaluate(originalData, evaluation);
           tsneCompressedData.setColumn(tsneCompressedData.getNcols()-1, evaluation);
           if ( config.getGeneralConfig().targetFileType == VisualizationFileType::CSV ) {
             CSVTools::writeMatrixToCSVFile(config.getGeneralConfig().currentDirectory +
               "/tsneCompression", tsneCompressedData);
           } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
             if (config.getVisualizationParameters().targetDimension != 2) {
               std::cout << "A json output is only available for compressions in 2 dimensions"
               "Storing the CSV instead" << std::endl;
               CSVTools::writeMatrixToCSVFile(config.getGeneralConfig().currentDirectory +
                 "/tsneCompression", tsneCompressedData);
             }
               storeTsneJson(tsneCompressedData, model);
             }
          }
        }
     }
    #pragma omp section
    {
      for(size_t index=0; index < models->size();index++) {

       if (config.getGeneralConfig().crossValidation) {
        config.getGeneralConfig().currentDirectory = config.getGeneralConfig().
                  targetDirectory+"/Fold_" + std::to_string(fold)
                  + "/Batch_" + std::to_string(batch)
                  +"/Model_" + std::to_string(index);
        } else {
         config.getGeneralConfig().currentDirectory = config.getGeneralConfig().
                   targetDirectory+"/Batch_" + std::to_string(batch)
                   +"/Model_" + std::to_string(index);
        }

        std::string mkdir("mkdir --parents ");

        mkdir.append(config.getGeneralConfig().currentDirectory);

        system(mkdir.data());

        auto currentModel = &(models->at(index));
        getLinearCuts(**currentModel);
        getHeatmap(**currentModel);
      }
    }
  }
}


void VisualizerClassification::getHeatmapsClassification(ModelFittingBase &model) {
  std::cout << "Generating the classification heatmaps" << std::endl;

  auto nDimensions = model.getDataset()->getDimension();

  if ( nDimensions == 1 ) {
    std::cout << "Heatmap generation is not available for models of 1 dimension" <<std::endl;
    return;
  }

  DataMatrix heatMapMatrix(0, nDimensions);
  if ( nDimensions >=3 ) {
    if ( nDimensions >= 4 ) {
     getHeatmapMore4DClassification(heatMapMatrix, model);
    } else if ( nDimensions == 3 ) {
      getHeatmap3D(heatMapMatrix, model);
    }
  } else {
   getHeatmap2DClassification(heatMapMatrix, model);
  }
}
void VisualizerClassification::getHeatmapMore4DClassification(DataMatrix &heatMapMatrix,
ModelFittingBase &model) {
  std::string outputDir(config.getGeneralConfig().currentDirectory+"/Classification"+
    "/Heatmaps_Classification/");

  std::cout << "Resolution " << std::to_string(resolution) << std::endl;
  double step = 1.0/resolution;

  for (double dim1 = 0; dim1 <= 1; dim1+=0.25) {
    for (double dim2 = 0; dim2 <= 1; dim2+=0.25) {
      for (double dim3 = 0; dim3 <= 1; dim3 = dim3+step) {
        for (double dim4 = 0; dim4 <= 1; dim4 = dim4+step) {
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

    std::string command("mkdir "+subfolder+" --parents");

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
           storeHeatmapJsonClassification(heatMapResults, model,
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
           storeHeatmapJsonClassification(heatMapResults, model,
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


void VisualizerClassification::getHeatmap2DClassification(DataMatrix &heatMapMatrix,
  ModelFittingBase &model) {
  std::string outputDir(config.getGeneralConfig().currentDirectory+"/");

  std::cout << "Resolution " << std::to_string(resolution) << std::endl;
  double step = 1.0/resolution;

  std::cout << "Step size: " << step<<std::endl;
  for (double dim1 = 0; dim1 <= 1; dim1 = dim1 + step) {
    for (double dim2 = 0; dim2 <= 1; dim2 = dim2 + step) {
      DataVector row(2, 0.5);
      row.set(0, dim1);
      row.set(1, dim2);
      heatMapMatrix.appendRow(row);
      // std::cout << row.toString() <<std::endl;
    }
  }

  DataMatrix heatMapResults(heatMapMatrix);
  DataVector evaluation(heatMapMatrix.getNrows());

  model.evaluate(heatMapMatrix, evaluation);

  heatMapResults.appendCol(evaluation);

  if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
  CSVTools::writeMatrixToCSVFile(outputDir+"ClassificationModel", heatMapResults);
  } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
   storeHeatmapJsonClassification(heatMapResults,
    model, outputDir+"ClassificationModel");
  }
}

void VisualizerClassification::storeTsneJson(DataMatrix &matrix, ModelFittingBase &model) {
  JSON jsonOutput;

  ModelFittingClassification* classificationModel =
    dynamic_cast<ModelFittingClassification*>(&model);

  auto models = classificationModel->getModels();

  size_t numberModels = models->size();

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

  DataVector classes;

  getClasses(zCol, models->size(), classes);

  jsonOutput["data"][0]["marker"].addIDAttr("color", zCol.toString());

  jsonOutput["data"][0]["marker"].addIDAttr("colorscale", "\"Viridis\"");

  jsonOutput["data"][0]["marker"].addIDAttr("opacity", 0.8);

  jsonOutput["data"][0]["marker"].addIDAttr("showscale", true);

  jsonOutput["data"][0]["marker"].addDictAttr("colorbar");

  jsonOutput["data"][0]["marker"]["colorbar"].addDictAttr("title");
  jsonOutput["data"][0]["marker"]["colorbar"]["title"].addIDAttr("text", "\"Class \"");
  jsonOutput["data"][0]["marker"]["colorbar"].addIDAttr("tickmode", "\"array\"");
  jsonOutput["data"][0]["marker"]["colorbar"].addIDAttr("tickvals", classes.toString());

  jsonOutput.addDictAttr("layout");

  jsonOutput["layout"].addDictAttr("title");

  jsonOutput["layout"]["title"].addIDAttr("text", "\"TSNE Compression\"");
  jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

  jsonOutput.serialize(config.getGeneralConfig().currentDirectory + "/tsneCompression.json");
}

void VisualizerClassification::storeHeatmapJsonClassification(DataMatrix &matrix, ModelFittingBase &model,
std::vector<size_t> indexes, size_t &varDim1, size_t &varDim2, std::string filepath) {
  indexes.erase(std::find(indexes.begin(), indexes.end(), varDim1));

  indexes.erase(std::find(indexes.begin(), indexes.end(), varDim2));

  ModelFittingClassification* classificationModel =
    dynamic_cast<ModelFittingClassification*>(&model);

  auto models = classificationModel->getModels();

  size_t numberModels = models->size();

  double maxValue = matrix.max(matrix.getNcols()-1);
  double minValue = matrix.min(matrix.getNcols()-1);

  JSON jsonOutput;
  jsonOutput.addListAttr("data");
  jsonOutput.addDictAttr("layout");
  jsonOutput["layout"].addDictAttr("title");

  if (matrix.getNcols() >= 4) {
    jsonOutput["layout"].addIDAttr("height", "1500");
  }
  jsonOutput["layout"].addListAttr("annotations");
  jsonOutput["layout"]["title"].addIDAttr("text", "\"Classification "
    "Heatmaps: Variable dimensions: " +
  std::to_string(varDim1 + 1) + " and " + std::to_string(varDim2 + 1) + "\"");
  jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

  DataVector zCol(matrix.getNrows());

  matrix.getColumn(matrix.getNcols()-1, zCol);

  DataVector classes;

  getClasses(zCol, models->size(), classes);

  size_t totalGraphs = ((matrix.getNcols()==4)?5:25);

  size_t rowsPerGraph = matrix.getNrows()/totalGraphs;

  unsigned int graphIndex = 0;

  for (unsigned int graphNumber = 0; graphNumber < totalGraphs; graphNumber++) {
    DataMatrix temp(matrix.getNrows(), matrix.getNcols());
    temp.copyFrom(matrix);

    unsigned int beginIndex = graphNumber*rowsPerGraph+1;

    temp.resizeToSubMatrix(beginIndex, 1, beginIndex+rowsPerGraph-1, matrix.getNcols());

    // Adding data trace
    jsonOutput["data"].addDictValue();

    jsonOutput["data"][graphIndex].addIDAttr("type", "\"heatmap\"");

    std::string xAxis("x" + ((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    jsonOutput["data"][graphIndex].addIDAttr("xaxis",
    "\"x"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");

    std::string yAxis("y" + ((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    jsonOutput["data"][graphIndex].addIDAttr("yaxis",
    "\"y"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");

    DataVector xCol(temp.getNrows());

    temp.getColumn(varDim1, xCol);

    jsonOutput["data"][graphIndex].addIDAttr("x", xCol.toString());

    DataVector yCol(temp.getNrows());

    temp.getColumn(varDim2, yCol);

    jsonOutput["data"][graphIndex].addIDAttr("y", yCol.toString());

    DataVector zCol(temp.getNrows());

    temp.getColumn(temp.getNcols()-1, zCol);

    jsonOutput["data"][graphIndex].addIDAttr("z", zCol.toString());
    jsonOutput["data"][graphIndex].addIDAttr("hoverinfo", "\"x+y+z\"");

    jsonOutput["data"][graphIndex].addIDAttr("zmin", minValue);
    jsonOutput["data"][graphIndex].addIDAttr("zmax", maxValue);
    jsonOutput["data"][graphIndex].addIDAttr("colorscale", "\"Viridis\"");

    if(graphNumber != 0 ) {
     jsonOutput["data"][graphIndex].addIDAttr("showscale", false);
    }
    else{
     jsonOutput["data"][graphIndex].addDictAttr("colorbar");
     jsonOutput["data"][graphIndex]["colorbar"].addIDAttr("tickmode", "\"array\"");
     jsonOutput["data"][graphIndex]["colorbar"].addIDAttr("tickvals", classes.toString());
    }

    graphIndex++;

    DataVector evaluation(originalData.getNrows());

    classificationModel->evaluate(originalData, evaluation);

    DataMatrix evaluatedData(originalData.getNrows(), originalData.getNcols());

    evaluatedData.copyFrom(originalData);

    evaluatedData.appendCol(evaluation);

    // Adding the data points
    for(size_t i = 0; i < numberModels ; i++) {
      DataMatrix tracePoints(0, evaluatedData.getNcols());
      for (size_t index = 0; index < evaluatedData.getNrows(); index++) {

        DataVector row(evaluatedData.getNcols());
        evaluatedData.getRow(index, row);
        if (row.get(row.getSize()-1) == classes.at(i)) {
          tracePoints.appendRow(row);
        }

      }
      jsonOutput["data"].addDictValue();
      jsonOutput["data"][graphIndex].addIDAttr("type", "\"scatter\"");
      jsonOutput["data"][graphIndex].addIDAttr("mode", "\"markers\"");


      DataVector xData(tracePoints.getNrows());

      tracePoints.getColumn(varDim1, xData);

      jsonOutput["data"][graphIndex].addIDAttr("x", xData.toString());

      DataVector yData(tracePoints.getNrows());

      tracePoints.getColumn(varDim2, yData);

      jsonOutput["data"][graphIndex].addIDAttr("y", yData.toString());
      jsonOutput["data"][graphIndex].addIDAttr("hoverinfo", "\"x+y+z\"");
      DataVector zData(tracePoints.getNrows());

      tracePoints.getColumn(tracePoints.getNcols()-1, zData);

      jsonOutput["data"][graphIndex].addDictAttr("marker");
      jsonOutput["data"][graphIndex]["marker"].
      addIDAttr("symbol", "\"star-dot\"");
      jsonOutput["data"][graphIndex]["marker"].
      addIDAttr("size", (long int)10);
      jsonOutput["data"][graphIndex]["marker"].
      addDictAttr("line");
      jsonOutput["data"][graphIndex]["marker"]["line"].
      addIDAttr("color", "\"blue\"");
      jsonOutput["data"][graphIndex]["marker"]["line"].
      addIDAttr("width", 1.5);
      jsonOutput["data"][graphIndex].
      addIDAttr("legendgroup", std::to_string(i));
      jsonOutput["data"][graphIndex]["marker"].
      addIDAttr("color", zData.toString());
      jsonOutput["data"][graphIndex]["marker"].
      addIDAttr("colorscale", "\"Viridis\"");
      jsonOutput["data"][graphIndex]["marker"].
      addIDAttr("cmin", matrix.min(matrix.getNcols()-1));
      jsonOutput["data"][graphIndex]["marker"].
      addIDAttr("cmax", matrix.max(matrix.getNcols()-1));
      jsonOutput["data"][graphIndex]["marker"].
      addIDAttr("showscale", false);

      jsonOutput["data"][graphIndex].addIDAttr("xaxis",
        "\"x"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");
      jsonOutput["data"][graphIndex].addIDAttr("yaxis",
        "\"y"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");

      if (graphNumber == 0) {
        jsonOutput["data"][graphIndex].addIDAttr("showlegend", true);
        jsonOutput["data"][graphIndex].addIDAttr("name", "\"Data Class  " +
          std::to_string((int)classes.at(i))+"\"");
      } else {
        jsonOutput["data"][graphIndex].addIDAttr("showlegend", false);
      }
      graphIndex++;
    }

    // Adding the grids
    DataVector firstRow(temp.getNcols());
    temp.getRow(0, firstRow);

    for(size_t i = 0; i < numberModels ; i++) {

      auto grid = models->at(i)->getGrid().clone();

      DataMatrix gridMatrix;

      grid->getStorage().getCoordinateArrays(gridMatrix);
      DataMatrix tempGrid(0, gridMatrix.getNcols());

      for (size_t index = 0; index < gridMatrix.getNrows(); index++) {
        DataVector row(gridMatrix.getNcols());
        gridMatrix.getRow(index, row);

        if (gridMatrix.getNcols() >= 4 ) {
          if ((row.get(indexes.at(0)) == firstRow.get(indexes.at(0)) &
          row.get(indexes.at(1)) == firstRow.get(indexes.at(1)))) {
            tempGrid.appendRow(row);
          }
        } else {
          if (row.get(indexes.at(0)) == firstRow.get(indexes.at(0))) {
            tempGrid.appendRow(row);
          }
        }
      }


      jsonOutput["data"].addDictValue();

      jsonOutput["data"][graphIndex].addIDAttr("type", "\"scatter\"");
      jsonOutput["data"][graphIndex].addIDAttr("mode", "\"markers\"");
      jsonOutput["data"][graphIndex].addIDAttr("xaxis",
        "\"x"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");
      jsonOutput["data"][graphIndex].addIDAttr("yaxis",
        "\"y"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");

      jsonOutput["data"][graphIndex].addIDAttr("legendgroup",
        std::to_string(i+numberModels));

      jsonOutput["data"][graphIndex].addDictAttr("marker");

      jsonOutput["data"][graphIndex]["marker"].
           addIDAttr("color", "\""+colors.at(i)+"\"");

      if (graphNumber == ((gridMatrix.getNcols()==3)?1:6)) {
        jsonOutput["data"][graphIndex].addIDAttr("showlegend", true);
        jsonOutput["data"][graphIndex].addIDAttr("name", "\" Grid Class " +
          std::to_string((int)classes.at(i))+"\"");
      } else {
        jsonOutput["data"][graphIndex].addIDAttr("showlegend", false);
      }

      DataVector xColGrid(tempGrid.getNrows());

      tempGrid.getColumn(varDim1, xColGrid);

      jsonOutput["data"][graphIndex].addIDAttr("x", xColGrid.toString());


      DataVector yColGrid(tempGrid.getNrows());

      tempGrid.getColumn(varDim2, yColGrid);
      jsonOutput["data"][graphIndex].addIDAttr("y", yColGrid.toString());
      jsonOutput["data"][graphIndex].addIDAttr("hoverinfo", "\"x+y\"");

      graphIndex++;

    }
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

    if (matrix.getNcols() > 4) {
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

    if (matrix.getNcols() >= 4) {
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

    if (matrix.getNcols() >= 4) {
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

    jsonOutput["layout"].addDictAttr("legend");
    jsonOutput["layout"]["legend"].addIDAttr("x", -0.2);
    jsonOutput["layout"]["legend"].addIDAttr("y", 1.0);


  std::cout << "Writing file " << filepath + ".json" << std::endl;
  jsonOutput.serialize(filepath+".json");
}


void VisualizerClassification::storeHeatmapJsonClassification(DataMatrix &matrix, ModelFittingBase &model,
std::string filepath) {


  ModelFittingClassification* classificationModel =
     dynamic_cast<ModelFittingClassification*>(&model);

  auto models = classificationModel->getModels();

  auto numberModels = models->size();

  unsigned int graphIndex = 0;

  JSON jsonOutput;
  jsonOutput.addListAttr("data");

  jsonOutput.addDictAttr("layout");

  jsonOutput["layout"].addDictAttr("title");

  jsonOutput["layout"].addListAttr("annotations");

  jsonOutput["layout"]["title"].addIDAttr("text", "\"Classification Heatmap\"");
  jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

  // Adding data trace
  jsonOutput["data"].addDictValue();

  jsonOutput["data"][graphIndex].addIDAttr("type", "\"heatmap\"");

  DataVector xCol(matrix.getNrows());

  matrix.getColumn(0, xCol);

  jsonOutput["data"][graphIndex].addIDAttr("x", xCol.toString());

  DataVector yCol(matrix.getNrows());

  matrix.getColumn(1, yCol);

  jsonOutput["data"][graphIndex].addIDAttr("y", yCol.toString());

  DataVector zCol(matrix.getNrows());

  matrix.getColumn(2, zCol);

  DataVector classes;

  getClasses(zCol, models->size(), classes);


  jsonOutput["data"][graphIndex].addIDAttr("z", zCol.toString());
  // jsonOutput["data"][0].addIDAttr("showlegend", false);
  jsonOutput["data"][graphIndex].addIDAttr("hoverinfo", "\"x+y+z\"");
  jsonOutput["data"][graphIndex].addIDAttr("colorscale", "\"Viridis\"");
  jsonOutput["data"][graphIndex].addIDAttr("connectgaps", true);
  jsonOutput["data"][graphIndex].addDictAttr("colorbar");
  jsonOutput["data"][graphIndex]["colorbar"].addIDAttr("tickmode", "\"array\"");
  jsonOutput["data"][graphIndex]["colorbar"].addIDAttr("tickvals", classes.toString());

  DataVector evaluation(originalData.getNrows());

  classificationModel->evaluate(originalData, evaluation);

  DataMatrix evaluatedData(originalData.getNrows(), originalData.getNcols());

  evaluatedData.copyFrom(originalData);

  evaluatedData.appendCol(evaluation);

  graphIndex++;
  // Adding the data points
  for(size_t i = 0; i < classes.size() ; i++) {

    DataMatrix tracePoints(0, evaluatedData.getNcols());
    for (size_t index = 0; index < evaluatedData.getNrows(); index++) {

      DataVector row(evaluatedData.getNcols());
      evaluatedData.getRow(index, row);
      if (row.get(row.getSize()-1) == classes.at(i)) {
        tracePoints.appendRow(row);
      }

    }
    jsonOutput["data"].addDictValue();
    jsonOutput["data"][graphIndex].addIDAttr("type", "\"scatter\"");
    jsonOutput["data"][graphIndex].addIDAttr("mode", "\"markers\"");


    DataVector xData(tracePoints.getNrows());

    tracePoints.getColumn(0, xData);

    jsonOutput["data"][graphIndex].addIDAttr("x", xData.toString());

    DataVector yData(tracePoints.getNrows());

    tracePoints.getColumn(1, yData);

    jsonOutput["data"][graphIndex].addIDAttr("y", yData.toString());
    jsonOutput["data"][graphIndex].addIDAttr("showlegend", true);
    jsonOutput["data"][graphIndex].addIDAttr("name", "\"Data Class " +
      std::to_string((int)classes.at(i))+"\"");
    jsonOutput["data"][graphIndex].addIDAttr("hoverinfo", "\"x+y+z\"");
    DataVector zData(tracePoints.getNrows());

    tracePoints.getColumn(2, zData);

    jsonOutput["data"][graphIndex].addDictAttr("marker");
    jsonOutput["data"][graphIndex]["marker"].addIDAttr("symbol", "\"star-dot\"");
    jsonOutput["data"][graphIndex]["marker"].addIDAttr("size", (long int)10);
    jsonOutput["data"][graphIndex]["marker"].addDictAttr("line");
    jsonOutput["data"][graphIndex]["marker"]["line"].addIDAttr("color", "\"blue\"");
    jsonOutput["data"][graphIndex]["marker"]["line"].addIDAttr("width", 1.5);
    jsonOutput["data"][graphIndex]["marker"].addIDAttr("color", zData.toString());
    jsonOutput["data"][graphIndex]["marker"].addIDAttr("colorscale", "\"Viridis\"");
    jsonOutput["data"][graphIndex]["marker"].addIDAttr("cmin", matrix.min(matrix.getNcols()-1));
    jsonOutput["data"][graphIndex]["marker"].addIDAttr("cmax", matrix.max(matrix.getNcols()-1));
    jsonOutput["data"][graphIndex]["marker"].addIDAttr("showscale", false);

    graphIndex++;
  }
  // Adding the sparse grids
  for(size_t i = 0; i < models->size() ; i++) {
    auto grid = models->at(i)->getGrid().clone();

    DataMatrix gridMatrix;

    grid->getStorage().getCoordinateArrays(gridMatrix);

    jsonOutput["data"].addDictValue();
    jsonOutput["data"][graphIndex].addIDAttr("type", "\"scatter\"");
    jsonOutput["data"][graphIndex].addIDAttr("mode", "\"markers\"");
    jsonOutput["data"][graphIndex].addDictAttr("marker");

    DataVector xColGrid(gridMatrix.getNrows());

    gridMatrix.getColumn(0, xColGrid);

    jsonOutput["data"][graphIndex].addIDAttr("x", xColGrid.toString());

    DataVector yColGrid(gridMatrix.getNrows());

    gridMatrix.getColumn(1, yColGrid);

    jsonOutput["data"][graphIndex].addIDAttr("y", yColGrid.toString());
    jsonOutput["data"][graphIndex].addIDAttr("showlegend", true);
    jsonOutput["data"][graphIndex].addIDAttr("name", "\"Grid Class " +
      std::to_string((int)classes.at(i))+"\"");
    jsonOutput["data"][graphIndex].addIDAttr("hoverinfo", "\"x+y\"");

    graphIndex++;
  }

  jsonOutput["layout"].addDictAttr("legend");
  jsonOutput["layout"]["legend"].addIDAttr("x", -0.15);
  jsonOutput["layout"]["legend"].addIDAttr("y", 1.0);

  std::cout << "Writing file " << filepath + ".json" << std::endl;
  jsonOutput.serialize(filepath + ".json");
}

void VisualizerClassification::
getClasses(DataVector &column, size_t numClasses, DataVector &classes) {
  size_t index = 0;
  while(classes.getSize() < numClasses) {
    int classValue = std::round(column.get(index));

     if (std::find(classes.begin(), classes.end(), classValue) == classes.end()) {
       classes.append(classValue);
     }
     index++;
  }
}

} // namespace datadriven
} // namespace sgpp
