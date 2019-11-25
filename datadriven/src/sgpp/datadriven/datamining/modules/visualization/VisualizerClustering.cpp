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


namespace sgpp {
namespace datadriven {
  VisualizerClustering::VisualizerClustering(VisualizerConfiguration config) {
    this->config = config;
  }

  void VisualizerClustering::runVisualization(ModelFittingBase &model, DataSource &dataSource,
      size_t fold, size_t batch) {
    model.getDataset()->getDimension();
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


    omp_set_num_threads(static_cast<int> (config.getVisualizationParameters().numberCores));


    #pragma omp parallel sections
    {
      #pragma omp section
      {
        std::string outputDirectory = currentDirectory + "/Clustering";
        createFolder(outputDirectory);

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
         // VisualizerClassification::getHeatmapsClassification(**classificationModel,
          //    outputDirectory, heatMapMatrixCls);

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

    }
  }
}  // namespace datadriven
}  // namespace sgpp
