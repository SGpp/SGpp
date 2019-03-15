/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseGridMinerCrossValidation.cpp
 *
 *  Created on: Jul 26, 2018
 *      Author: dominik
 */

// FOR EVALUATION ONLY
#include <sgpp/datadriven/datamining/base/evaluationtools.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorFactory.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerCrossValidation.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <iostream>
#include <vector>
using evalu::getTime;

namespace sgpp {
namespace datadriven {

SparseGridMinerCrossValidation::SparseGridMinerCrossValidation(
    DataSourceCrossValidation* dataSource, ModelFittingBase* fitter, Scorer* scorer)
    : SparseGridMiner(fitter, scorer), dataSource{dataSource} {}

double SparseGridMinerCrossValidation::learn(bool verbose) {
  // todo(fuchsgdk): see below

  const CrossvalidationConfiguration& crossValidationConfig =
      dataSource->getCrossValidationConfig();

  std::vector<double> scores;
  scores.reserve(crossValidationConfig.kfold_);

  for (size_t fold = 0; fold < crossValidationConfig.kfold_; fold++) {
    dataSource->setFold(fold);

    // todo(fuchsgdk):
    // This is the kind of cv implemented by Lettrich in the scorer class and it was
    // merely moved to fit into the data source. Conceptual changes might be done in order to
    // really support batch based learning with cv and not only regression.
    // What should be done is reimplementing the data source such that it provides batches

    std::cout << "###############"
              << "Fold #" << fold << std::endl;

    // Create a refinement monitor for this fold
    RefinementMonitorFactory monitorFactory;
    RefinementMonitor* monitor = monitorFactory.createRefinementMonitor(
        fitter->getFitterConfiguration().getRefinementConfig());

    // Reset the fitter
    fitter->reset();

    for (size_t epoch = 0; epoch < dataSource->getConfig().epochs; epoch++) {
      std::cout << "###############"
                << "Starting training epoch #" << epoch << std::endl;
      dataSource->reset();
      Dataset* validationData = dataSource->getValidationData();
      size_t validationSize = validationData->getNumberInstances();
      std::cout << "Validation data size: " << validationSize << std::endl;
      // Process dataset iteratively
      size_t iteration = 0;
      while (true) {
        std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
        size_t numInstances = dataset->getNumberInstances();
        if (numInstances == 0) {
          // The source does not provide any more samples
          break;
        }
        std::cout << "###############" << getTime() << "Itertation #" << (iteration++) << std::endl
                  << "Batch size: " << numInstances << std::endl;

        // Train model on new batch
        fitter->update(*dataset);

        // Evaluate the score on the training and validation data
        double scoreTrain = scorer->test(*fitter, *dataset);
        double scoreVal = scorer->test(*fitter, *validationData);

        std::cout << getTime() << std::endl
                  << "Score on batch: " << scoreTrain << std::endl
                  << "Score on validation data: " << scoreVal << std::endl;

        // Refine the model if neccessary
        monitor->pushToBuffer(numInstances, scoreVal, scoreTrain);
        size_t refinements = monitor->refinementsNecessary();
        while (refinements--) {
          std::cout << getTime();
          fitter->refine();
          std::cout << getTime();
        }

        std::cout << "###############" << getTime() << "Iteration finished." << std::endl;
      }
    }
    // Evaluate the final score on the validation data
    dataSource->reset();
    Dataset* validationData = dataSource->getValidationData();
    scores.push_back(scorer->test(*fitter, *validationData));
  }

  // Calculate mean score and std deviation
  double meanScore = 0.0;
  for (size_t idx = 0; idx < scores.size(); idx++) {
    meanScore += scores[idx];
  }
  meanScore /= static_cast<double>(scores.size());
  double stdDeviation = 0.0;
  for (size_t idx = 0; idx < scores.size(); idx++) {
    stdDeviation += std::pow(scores[idx] - meanScore, 2);
  }
  stdDeviation = std::sqrt(stdDeviation / static_cast<double>(crossValidationConfig.kfold_ - 1));
  std::cout << "###############" << std::endl
            << "Mean score: " << meanScore << std::endl
            << "Standard deviation: " << stdDeviation << std::endl;
  return meanScore;
}
} /* namespace datadriven */
} /* namespace sgpp */
