// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/datadriven/datamining/modules/visualization/VisualizerClassification.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/ListNode.hpp>
#include <sgpp/base/tools/json/DictNode.hpp>
#include <omp.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>


namespace sgpp {
namespace datadriven {
VisualizerClassification::VisualizerClassification(VisualizerConfiguration config) {
  this->config = config;
}


void VisualizerClassification::runVisualization(ModelFittingBase &model, DataSource &dataSource,
  size_t fold, size_t batch) {
  if (batch % config.getGeneralConfig().numBatches_ != 0 ||
    !config.getGeneralConfig().execute_) {
    return;
  }

  if (fold == 0 && batch == 0) {
    originalData = dataSource.getAllSamples()->getData();
    resolution = static_cast<size_t>(pow(2,
      model.getFitterConfiguration().getGridConfig().level_+2));
  }

  ModelFittingClassification* classificationModel =
    dynamic_cast<ModelFittingClassification*>(&model);

  auto models = classificationModel->getModels();
  auto classIdx = classificationModel->getClassIdx();

  classes.resizeZero(models->size());

  for (auto const& x : classIdx) {
    classes.set(x.second, x.first);
  }

  createFolder(config.getGeneralConfig().
    targetDirectory_);
  // Creating the output directory
  if (config.getGeneralConfig().crossValidation_) {
    currentDirectory = config.getGeneralConfig().
    targetDirectory_+"/Fold_" + std::to_string(fold);
    createFolder(currentDirectory);
    currentDirectory = config.getGeneralConfig().
        targetDirectory_+"/Fold_" + std::to_string(fold) + "/Batch_" + std::to_string(batch);
    createFolder(currentDirectory);

  } else {
    currentDirectory = config.getGeneralConfig().
    targetDirectory_+"/Batch_" + std::to_string(batch);
    createFolder(currentDirectory);
  }

  std::cout << "Creating output directory " << config.getGeneralConfig().targetDirectory_
     << std::endl;

  createFolder(currentDirectory+"/Classification");

  omp_set_num_threads(static_cast<int> (config.getVisualizationParameters().numberCores_));

  #pragma omp parallel sections
  {
    #pragma omp section
    {
      if (std::find(config.getGeneralConfig().algorithm_.begin(),
                   config.getGeneralConfig().algorithm_.end(), "heatmaps") !=
                   config.getGeneralConfig().algorithm_.end()) {
        DataMatrix heatMapClassificationMatrix;
        initializeMatrices(model, heatMapClassificationMatrix);
        getHeatmapsClassification(model, currentDirectory, heatMapClassificationMatrix);
      }
    }
    /* To be added after tsne is functioning again
     * #pragma omp section
    {
      if (config.getGeneralConfig().algorithm == "tsne") {
        if (fold == 0 && batch == 0) {
          runTsne(model);
        }
        if (originalData.getNcols() >= 1) {
          DataVector evaluation(originalData.getNrows());
          model.evaluate(originalData, evaluation);
          tsneCompressedData.setColumn(tsneCompressedData.getNcols()-1, evaluation);
          if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
            CSVTools::writeMatrixToCSVFile(currentDirectory +
              "/tsneCompression", tsneCompressedData);
          } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
            if (config.getVisualizationParameters().targetDimension != 2) {
              std::cout << "A json output is only available for compressions in 2 dimensions"
              "Storing the CSV instead" << std::endl;
              CSVTools::writeMatrixToCSVFile(currentDirectory +
                "/tsneCompression", tsneCompressedData);
            }
              storeTsneJson(tsneCompressedData, model, currentDirectory);
          }
        }
      }
    }*/
    #pragma omp section
    {
      // Running the density estimation visualization for each model
      #pragma omp parallel for schedule(dynamic)
      for (size_t index=0; index < models->size(); index++) {
        DataMatrix heatMapMatrixThread;
        DataMatrix cutMatrixThread;
        VisualizerDensityEstimation::initializeMatrices(model, cutMatrixThread,
          heatMapMatrixThread);
        if (config.getGeneralConfig().crossValidation_) {
          currentDirectory = config.getGeneralConfig().
                   targetDirectory_+"/Fold_" + std::to_string(fold)
                   + "/Batch_" + std::to_string(batch)
                   +"/Model_Class_" + std::to_string(static_cast<int>(classes.get(index)));
         } else {
           currentDirectory = config.getGeneralConfig().
                    targetDirectory_+"/Batch_" + std::to_string(batch)
                    +"/Model_Class_" + std::to_string(static_cast<int>(classes.get(index)));
         }

        createFolder(currentDirectory);

        auto currentModel = &(models->at(index));
        #pragma omp parallel sections
        {
          #pragma omp section
          {
            if (std::find(config.getGeneralConfig().algorithm_.begin(),
                         config.getGeneralConfig().algorithm_.end(), "linearcuts")
                         != config.getGeneralConfig().algorithm_.end()) {
              getLinearCuts(**currentModel, currentDirectory, cutMatrixThread);
            }
          }
          #pragma omp section
          {
            if (std::find(config.getGeneralConfig().algorithm_.begin(),
                         config.getGeneralConfig().algorithm_.end(), "heatmaps")
                         != config.getGeneralConfig().algorithm_.end()) {
              getHeatmap(**currentModel, currentDirectory, heatMapMatrixThread);
            }
          }
        }
      }
    }
  }
}


void VisualizerClassification::initializeMatrices(ModelFittingBase &model,
  DataMatrix &classMatrix) {
    std::cout << "Resolution " << std::to_string(resolution) << std::endl;
    double step = 1.0/static_cast<double>(resolution);

    auto nDimensions = model.getDataset()->getDimension();

    classMatrix.resize(0, nDimensions);
    if (nDimensions >=3) {
      if (nDimensions >= 4) {
        for (double dim1 = 0; dim1 <= 1; dim1+=0.25) {
          for (double dim2 = 0; dim2 <= 1; dim2+=0.25) {
            for (double dim3 = 0; dim3 <= 1; dim3+= step) {
              for (double dim4 = 0; dim4 <= 1; dim4+= step) {
              DataVector row(classMatrix.getNcols(), 0.5);
              row.set(0, dim3);
              row.set(1, dim4);
              row.set(2, dim1);
              row.set(3, dim2);

              classMatrix.appendRow(row);
              }
            }
          }
        }
     } else {
       for (double dim1 = 0; dim1 <= 1; dim1+=0.25) {
         for (double dim2 = 0; dim2 <= 1; dim2 += step) {
           for (double dim3 = 0; dim3 <= 1; dim3 += step) {
             DataVector row(3, 0.5);
             row.set(0, dim3);
             row.set(1, dim2);
             row.set(2, dim1);
             classMatrix.appendRow(row);
           }
         }
       }
     }
    } else {
      for (double dim1 = 0; dim1 <= 1; dim1 += step) {
        for (double dim2 = 0; dim2 <= 1; dim2 += step) {
          DataVector row(2, 0.5);
          row.set(0, dim1);
          row.set(1, dim2);
          classMatrix.appendRow(row);
        }
      }
    }
}

void VisualizerClassification::getHeatmapsClassification(ModelFittingBase &model,
  std::string currentDirectory, DataMatrix &classMatrix) {
  std::cout << "Generating the classification heatmaps" << std::endl;

  auto nDimensions = model.getDataset()->getDimension();

  if ( nDimensions == 1 ) {
    std::cout << "Heatmap generation is not available for models of 1 dimension" <<std::endl;
    return;
  }
  if ( nDimensions >=3 ) {
    std::string outputDir(currentDirectory+"/Classification"+
       "/Heatmaps/");
    createFolder(outputDir);
    if ( nDimensions >= 4 ) {
     getHeatmapMore4DClassification(model, currentDirectory, classMatrix);
    } else {
      getHeatmap3DClassification(model, currentDirectory, classMatrix);
    }
  } else {
    getHeatmap2DClassification(model, currentDirectory, classMatrix);
  }
}

// The algorithm used for the heatmaps is this:
// 1° Evaluate the initial matrix
// 2° Store current evaluation
// 3° Shift to the right based on the column indexes
// 4° Repeat steps 1 to 3 two more times
// 5° Repeat the steps 1 to 4 but shifting to the left instead
// 6° Shift to the right if its the first time you reach this step
// 7° Shift to the left
// 8° Update the indexes
// 9° Repeat until the last column index exceeds the number of dimensions
void VisualizerClassification::getHeatmapMore4DClassification(
ModelFittingBase &model, std::string currentDirectory, DataMatrix &classMatrix) {
  std::string outputDir(currentDirectory+"/Classification"+
    "/Heatmaps/");

  std::vector <size_t> variableColumnIndexes = {0, 1, 2, 3};

  std::vector<size_t> workingIndexes = {0, 0, 0};

  while (variableColumnIndexes.at(3) < classMatrix.getNcols()) {
    std::string subfolder(outputDir+"dimensions_"+
    std::to_string(variableColumnIndexes.at(0)+1)+"_"+
    std::to_string(variableColumnIndexes.at(1)+1)+"_"+
    std::to_string(variableColumnIndexes.at(2)+1)+"_"+
    std::to_string(variableColumnIndexes.at(3)+1));

    createFolder(subfolder);
    std::copy(variableColumnIndexes.begin()+1, variableColumnIndexes.end(),
    workingIndexes.begin());

    for (size_t iteration = 0; iteration < 2; iteration++) {
      for (size_t combination = 0; combination < 3; combination++) {
        DataMatrix heatMapResults(classMatrix);
        DataVector evaluation(classMatrix.getNrows());

        model.evaluate(classMatrix, evaluation);

        heatMapResults.appendCol(evaluation);

        if (iteration == 0) {
          if (config.getGeneralConfig().targetFileType_ == VisualizationFileType::CSV) {
            CSVTools::writeMatrixToCSVFile(subfolder+"/Heatmap_var_dimensions_"
            +std::to_string(variableColumnIndexes.at(0)+1)+"_"+
            std::to_string(variableColumnIndexes.at(combination+1)+1), heatMapResults);
          } else if (config.getGeneralConfig().targetFileType_ == VisualizationFileType::json) {
           storeHeatmapJsonClassification(heatMapResults, model,
            variableColumnIndexes, variableColumnIndexes.at(0),
            variableColumnIndexes.at(combination+1),
            subfolder+"/Heatmap_var_dimensions_"
                       +std::to_string(variableColumnIndexes.at(0)+1)+"_"+
                       std::to_string(variableColumnIndexes.at(combination+1)+1));
          }
          translateColumnsRight(classMatrix, workingIndexes);
        } else {
          if (config.getGeneralConfig().targetFileType_ == VisualizationFileType::CSV) {
            CSVTools::writeMatrixToCSVFile(subfolder+"/Heatmap_var_dimensions_"
              +std::to_string(variableColumnIndexes.at((combination < 2)?1:2)+1)+"_"+
              std::to_string(variableColumnIndexes.at((combination < 1)?2:3)+1), heatMapResults);
          } else if (config.getGeneralConfig().targetFileType_ == VisualizationFileType::json) {
           storeHeatmapJsonClassification(heatMapResults, model,
             variableColumnIndexes, variableColumnIndexes.at((combination < 2)?1:2),
             variableColumnIndexes.at((combination < 1)?2:3),
             subfolder+"/Heatmap_var_dimensions_"
             +std::to_string(variableColumnIndexes.at((combination < 2)?1:2)+1)+"_"+
             std::to_string(variableColumnIndexes.at((combination < 1)?2:3)+1));
          }
          translateColumnsLeft(classMatrix, workingIndexes);
        }
      }

      if (iteration == 0) {
        translateColumnsRight(classMatrix, variableColumnIndexes);
      }
    }
    translateColumnsLeft(classMatrix, variableColumnIndexes);

    updateIndexesHeatmap(variableColumnIndexes, classMatrix);
  }
}

// The algorithm used for the heatmaps is this:
// 1° Evaluate the initial matrix
// 2° Store current evaluation
// 3° Shift to the right
// 4° Repeat steps 1 to 3 two more times
void VisualizerClassification::getHeatmap3DClassification(ModelFittingBase &model,
  std::string currentDirectory, DataMatrix &classMatrix) {
  std::string outputDir(currentDirectory+"/Classification/");

  // Dummy to reutilize the storejson method
  std::vector <size_t> variableColumnIndexes = {0, 1, 2};

  for (size_t combination = 0; combination < 3; combination++) {
    DataMatrix heatMapResults(classMatrix);
    DataVector evaluation(classMatrix.getNrows());

    model.evaluate(classMatrix, evaluation);

    heatMapResults.appendCol(evaluation);

    translateColumns(classMatrix, classMatrix.getNcols());
    if (config.getGeneralConfig().targetFileType_ == VisualizationFileType::CSV) {
      CSVTools::writeMatrixToCSVFile(outputDir+"Heatmap_var_dimensions_"
      +std::to_string(combination+1)
      +"_"+((combination < 2)?std::to_string(combination+2):std::to_string(1)), heatMapResults);

    } else if (config.getGeneralConfig().targetFileType_ == VisualizationFileType::json) {
      storeHeatmapJsonClassification(heatMapResults,
      model,
      variableColumnIndexes,
      variableColumnIndexes.at(combination),
      variableColumnIndexes.at((combination < 2)?combination+1:0),
      outputDir+"Heatmap_var_dimensions_"+std::to_string(combination+1)+"_"+
      ((combination < 2)?std::to_string(combination+2):std::to_string(1)));
    }
  }
}
void VisualizerClassification::getHeatmap2DClassification(
  ModelFittingBase &model, std::string currentDirectory, DataMatrix &classMatrix) {
  std::string outputDir(currentDirectory+"/Classification/");

  DataMatrix heatMapResults(classMatrix);
  DataVector evaluation(classMatrix.getNrows());

  model.evaluate(classMatrix, evaluation);

  heatMapResults.appendCol(evaluation);

  if (config.getGeneralConfig().targetFileType_ == VisualizationFileType::CSV) {
  CSVTools::writeMatrixToCSVFile(outputDir+"ClassificationModel", heatMapResults);
  } else if (config.getGeneralConfig().targetFileType_ == VisualizationFileType::json) {
    storeHeatmapJsonClassification(heatMapResults,
    model, outputDir+"ClassificationModel");
  }
}

void VisualizerClassification::storeTsneJson(DataMatrix &matrix, ModelFittingBase &model,
  std::string currentDirectory) {
  json::JSON jsonOutput;

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
  jsonOutput["data"][0]["marker"]["colorbar"]["title"].addIDAttr("text", "\"Class \"");
  jsonOutput["data"][0]["marker"]["colorbar"].addIDAttr("tickmode", "\"array\"");
  jsonOutput["data"][0]["marker"]["colorbar"].addIDAttr("tickvals", classes.toString());

  jsonOutput.addDictAttr("layout");

  jsonOutput["layout"].addDictAttr("title");

  jsonOutput["layout"]["title"].addIDAttr("text", "\"TSNE Compression\"");
  jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

  jsonOutput.serialize(currentDirectory + "/Classification/tsneCompression.json");
}

void VisualizerClassification::storeHeatmapJsonClassification(DataMatrix &matrix,
  ModelFittingBase &model, std::vector<size_t> indexes, size_t &varDim1,
  size_t &varDim2, std::string filepath) {
  indexes.erase(std::find(indexes.begin(), indexes.end(), varDim1));

  indexes.erase(std::find(indexes.begin(), indexes.end(), varDim2));

  ModelFittingClassification* classificationModel =
    dynamic_cast<ModelFittingClassification*>(&model);

  auto models = classificationModel->getModels();

  size_t numberModels = models->size();

  std::vector<bool> showLegendGroup(numberModels, true);

  double maxValue = matrix.max(matrix.getNcols()-1);
  double minValue = matrix.min(matrix.getNcols()-1);

  json::JSON jsonOutput;
  jsonOutput.addListAttr("data");
  jsonOutput.addDictAttr("layout");
  jsonOutput["layout"].addDictAttr("title");

  if (matrix.getNcols() > 4) {
    jsonOutput["layout"].addIDAttr("height", "1500");
  }
  jsonOutput["layout"].addListAttr("annotations");
  jsonOutput["layout"]["title"].addIDAttr("text", "\"Classification "
    "Heatmaps: Variable dimensions: " +
  std::to_string(varDim1 + 1) + " and " + std::to_string(varDim2 + 1) + "\"");
  jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

  DataVector zCol(matrix.getNrows());

  matrix.getColumn(matrix.getNcols()-1, zCol);

  size_t totalGraphs = ((matrix.getNcols() == 4)?5:25);

  size_t rowsPerGraph = matrix.getNrows()/totalGraphs;

  size_t graphIndex = 0;

  for (size_t graphNumber = 0; graphNumber < totalGraphs; graphNumber++) {
    DataMatrix temp(matrix.getNrows(), matrix.getNcols());
    temp.copyFrom(matrix);

    size_t beginIndex = graphNumber*rowsPerGraph+1;

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

    if (graphNumber != 0) {
      jsonOutput["data"][graphIndex].addIDAttr("showscale", false);
    } else {
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
    for (size_t i = 0; i < numberModels ; i++) {
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
      addIDAttr("size", static_cast<size_t>(10));
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
          std::to_string(static_cast<int>(classes.at(i)))+"\"");
      } else {
        jsonOutput["data"][graphIndex].addIDAttr("showlegend", false);
      }
      graphIndex++;
    }

    // Adding the grids
    DataVector firstRow(temp.getNcols());
    temp.getRow(0, firstRow);

    for (size_t i = 0; i < numberModels ; i++) {
      auto grid = models->at(i)->getGrid().clone();

      DataMatrix gridMatrix;

      grid->getStorage().getCoordinateArrays(gridMatrix);
      DataMatrix tempGrid(0, gridMatrix.getNcols());

      for (size_t index = 0; index < gridMatrix.getNrows(); index++) {
        DataVector row(gridMatrix.getNcols());
        gridMatrix.getRow(index, row);

        if (gridMatrix.getNcols() >= 4) {
          if (((row.get(indexes.at(0)) == firstRow.get(indexes.at(0))) &
          (row.get(indexes.at(1)) == firstRow.get(indexes.at(1))))) {
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

      if (showLegendGroup.at(i) && tempGrid.getNrows() > 0) {
        jsonOutput["data"][graphIndex].addIDAttr("showlegend", true);
        jsonOutput["data"][graphIndex].addIDAttr("name", "\" Grid Class " +
          std::to_string(static_cast<int>(classes.at(i)))+"\"");
        showLegendGroup.at(i) = false;
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
    jsonOutput["layout"][xAxisName]["domain"].addIdValue(0.2*
      static_cast<double>(graphNumber%5));
    jsonOutput["layout"][xAxisName]["domain"].addIdValue(0.2*
      (static_cast<double>(graphNumber%5)+1)-0.05);


    jsonOutput["layout"].addDictAttr(yAxisName);
    jsonOutput["layout"][yAxisName].addIDAttr("anchor", "\""+xAxis+"\"");
    jsonOutput["layout"][yAxisName].addIDAttr("type", "\"linear\"");

    if (matrix.getNcols() > 4) {
      jsonOutput["layout"][yAxisName].addListAttr("domain");
      jsonOutput["layout"][yAxisName]["domain"].addIdValue(0.2*
        static_cast<double>(graphNumber/5));
      jsonOutput["layout"][yAxisName]["domain"].addIdValue(0.2*
        (static_cast<double>(graphNumber/5)+1)-0.1);
    }

    // Adding titles to subplots
    std::string dim1Text(std::to_string(indexes.at(0)+1));
    std::string dim1ValueText(std::to_string(firstRow.get(indexes.at(0))));
    dim1ValueText.erase(dim1ValueText.find_last_not_of('0') + 2, std::string::npos);
    dim1Text = "\"Dim " +
    dim1Text + "=" + dim1ValueText;

    std::string dim2Text = "";
    std::string dim2ValueText = "\"";

    if (matrix.getNcols() > 4) {
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


    if (matrix.getNcols() > 4) {
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


void VisualizerClassification::storeHeatmapJsonClassification(DataMatrix &matrix,
  ModelFittingBase &model, std::string filepath) {
  ModelFittingClassification* classificationModel =
     dynamic_cast<ModelFittingClassification*>(&model);

  auto models = classificationModel->getModels();

  unsigned int graphIndex = 0;

  json::JSON jsonOutput;
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
  for (size_t i = 0; i < classes.size() ; i++) {
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
      std::to_string(static_cast<int>(classes.at(i)))+"\"");
    jsonOutput["data"][graphIndex].addIDAttr("hoverinfo", "\"x+y+z\"");
    DataVector zData(tracePoints.getNrows());

    tracePoints.getColumn(2, zData);

    jsonOutput["data"][graphIndex].addDictAttr("marker");
    jsonOutput["data"][graphIndex]["marker"].addIDAttr("symbol", "\"star-dot\"");
    jsonOutput["data"][graphIndex]["marker"].addIDAttr("size",
      static_cast<size_t>(10));
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
  for (size_t i = 0; i < models->size(); i++) {
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
      std::to_string(static_cast<int>(classes.at(i)))+"\"");
    jsonOutput["data"][graphIndex].addIDAttr("hoverinfo", "\"x+y\"");

    graphIndex++;
  }

  jsonOutput["layout"].addDictAttr("legend");
  jsonOutput["layout"]["legend"].addIDAttr("x", -0.15);
  jsonOutput["layout"]["legend"].addIDAttr("y", 1.0);

  std::cout << "Writing file " << filepath + ".json" << std::endl;
  jsonOutput.serialize(filepath + ".json");
}

}  // namespace datadriven
}  // namespace sgpp
