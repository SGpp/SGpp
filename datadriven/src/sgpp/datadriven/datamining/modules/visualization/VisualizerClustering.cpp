// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/datadriven/datamining/modules/visualization/VisualizerClustering.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClustering.hpp>

#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/ListNode.hpp>
#include <sgpp/base/tools/json/DictNode.hpp>

#include <omp.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <queue>


namespace sgpp {
namespace datadriven {
  VisualizerClustering::VisualizerClustering(VisualizerConfiguration config) {
    this->config = config;
  }

  void VisualizerClustering::runVisualization(ModelFittingBase &model, DataSource &dataSource,
      size_t fold, size_t batch) {
    /*model.getDataset()->getDimension();
    if (batch % config.getGeneralConfig().numBatches != 0 ||
        !config.getGeneralConfig().execute) {
      return;
    }

    if (fold == 0 && batch == 0) {
      originalData = dataSource.getAllSamples()->getData();
      resolution = static_cast<size_t>(pow(2,
                                           model.getFitterConfiguration().getGridConfig().level_ +
                                           2));
    }

    ModelFittingClustering *clusteringModel =
        dynamic_cast<ModelFittingClustering *>(&model);

    std::cout << "Creating output directory " << config.getGeneralConfig().targetDirectory
              << std::endl;

    createFolder(config.getGeneralConfig().
        targetDirectory);
    // Creating the output directory
    if (config.getGeneralConfig().crossValidation) {
      currentDirectory = config.getGeneralConfig().
          targetDirectory + "/Fold_" + std::to_string(fold);
      createFolder(currentDirectory);
      currentDirectory = config.getGeneralConfig().
          targetDirectory + "/Fold_" + std::to_string(fold) + "/Batch_" + std::to_string(batch);
      createFolder(currentDirectory);

    } else {
      currentDirectory = config.getGeneralConfig().
          targetDirectory + "/Batch_" + std::to_string(batch);
      createFolder(currentDirectory);
    }

    std::string outputDirectory = currentDirectory + "/Clustering";
    createFolder(outputDirectory);

    omp_set_num_threads(static_cast<int> (config.getVisualizationParameters().numberCores));

    DataMatrix points = clusteringModel->getLabeledPoints();
    storeTsneJson(points, model, outputDirectory);
    #pragma omp parallel sections
    {
      #pragma omp section
      {
        auto densityEstimationModel = clusteringModel->getDensityEstimationModel();
        auto classificationModel = clusteringModel->getClassificationModel();

        DataMatrix heatMapMatrixDE;
        DataMatrix cutMatrixDE;
        VisualizerDensityEstimation::initializeMatrices(model, cutMatrixDE,
                                                        heatMapMatrixDE);

        DataMatrix heatMapMatrixCls;
        VisualizerClassification::initializeMatrices(model, heatMapMatrixCls);

        auto models = (*classificationModel)->getModels();
        auto classIdx = (*classificationModel)->getClassIdx();

        classes.resizeZero(models->size());

        for (auto const& x : classIdx) {
          classes.set(x.second, x.first);
        }

        if (std::find(config.getGeneralConfig().algorithm.begin(),
                      config.getGeneralConfig().algorithm.end(), "heatmaps") !=
            config.getGeneralConfig().algorithm.end()) {
          VisualizerDensityEstimation::getHeatmap(**densityEstimationModel, outputDirectory,
                                                  heatMapMatrixDE);
          VisualizerClassification::getHeatmapsClassification(**classificationModel,
                outputDirectory, heatMapMatrixCls);
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
      }
    } */
  }

  void VisualizerClustering::storeTsneJson(DataMatrix &matrix, ModelFittingBase &model,
      std::string currentDirectory) {
    json::JSON jsonOutput;

    ModelFittingClustering *clusteringModel =
        dynamic_cast<ModelFittingClustering *>(&model);

    jsonOutput.addListAttr("data");

    // Trace for the data points
    jsonOutput["data"].addDictValue();
    jsonOutput["data"][0].addIDAttr("type", "\"scatter\"");
    jsonOutput["data"][0].addIDAttr("mode", "\"markers\"");

    DataVector xCol(matrix.getNrows());

    matrix.getColumn(0, xCol);

    jsonOutput["data"][0].addIDAttr("x", xCol.toString());

    DataVector yCol(matrix.getNrows());

    matrix.getColumn(1, yCol);
    jsonOutput["data"][0].addIDAttr("y", yCol.toString());

    // Trace for the edges
    jsonOutput["data"].addDictValue();
    jsonOutput["data"][1].addIDAttr("type", "\"scatter\"");
    jsonOutput["data"][1].addIDAttr("mode", "\"lines\"");
    jsonOutput["data"][1].addDictAttr("line");
    jsonOutput["data"][1]["line"].addIDAttr("width", 0.5);


    DataVector source(matrix.getNcols());
    DataVector sink(matrix.getNcols());

    jsonOutput["data"][1].addListAttr("x");
    jsonOutput["data"][1].addListAttr("y");

    /*for (size_t index = 0; index < matrix.getNrows(); index++) {
      matrix.getRow(index, source);
      std::priority_queue<VpHeapItem> heap2(nearestNeighbors[index]);
      while (!heap2.empty()) {
        jsonOutput["data"][1]["x"].addIdValue(source.get(0));
        jsonOutput["data"][1]["y"].addIdValue(source.get(1));
        matrix.getRow(heap2.top().index, sink);
        jsonOutput["data"][1]["x"].addIdValue(sink.get(0));
        jsonOutput["data"][1]["y"].addIdValue(sink.get(1));
        heap2.pop();
        jsonOutput["data"][1]["x"].addIdValue("\"None\"\n");
        jsonOutput["data"][1]["y"].addIdValue("\"None\"\n");
      }
    }*/
    /*
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
    jsonOutput["data"][0]["marker"]["colorbar"].addIDAttr("tickvals", classes.toString());*/

    jsonOutput.addDictAttr("layout");

    jsonOutput["layout"].addDictAttr("title");

    jsonOutput["layout"]["title"].addIDAttr("text", "\"Nearest Neighbors Graph\"");
    jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

    jsonOutput.serialize(currentDirectory + "/Graph.json");
  }
}  // namespace datadriven
}  // namespace sgpp
