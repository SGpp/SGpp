// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/base/SparseGridMinerCrossValidation.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorFactory.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <iostream>
#include <vector>

namespace sgpp {
namespace datadriven {

SparseGridMinerCrossValidation::SparseGridMinerCrossValidation(
    DataSourceCrossValidation* dataSource, ModelFittingBase* fitter, Scorer* scorer,
    Visualizer* visualizer)
    : SparseGridMiner(fitter, scorer, visualizer), dataSource{dataSource} {}

double SparseGridMinerCrossValidation::learn(bool verbose) {
// todo(fuchsgdk): see below

#ifdef USE_SCALAPACK
  if (fitter->getFitterConfiguration().getParallelConfig().scalapackEnabled_) {
    auto processGrid = fitter->getProcessGrid();
    if (!processGrid->isProcessInGrid()) {
      return 0.0;
    }
  }
#endif /* USE_SCALAPACK */

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

    std::ostringstream out;
    out << "###############"
        << "Fold #" << fold;
    print(out);

    // Create a refinement monitor for this fold
    RefinementMonitorFactory monitorFactory;
    RefinementMonitor* monitor = monitorFactory.createRefinementMonitor(
        fitter->getFitterConfiguration().getRefinementConfig());

    // Reset the fitter
    fitter->reset();

    for (size_t epoch = 0; epoch < dataSource->getConfig().epochs_; epoch++) {
      if (verbose) {
        std::ostringstream out;
        out << "###############"
            << "Starting training epoch #" << epoch;
        print(out);
      }
      dataSource->reset();
      Dataset* validationData = dataSource->getValidationData();
      size_t validationSize = validationData->getNumberInstances();

      if (verbose) {
        std::ostringstream out;
        out << "Validation data size: " << validationSize;
        print(out);
      }
      // Process dataset iteratively
      size_t iteration = 0;
      while (true) {
        std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
        size_t numInstances = dataset->getNumberInstances();
        if (numInstances == 0) {
          // The source does not provide any more samples
          break;
        }

        if (verbose) {
          std::ostringstream out;
          out << "###############"
              << "Iteration #" << (iteration) << std::endl
              << "Batch size: " << numInstances;
          print(out);
        }

        // Train model on new batch
        fitter->update(*dataset);

        // Evaluate the score on the training and validation data
        double scoreTrain = scorer->test(*fitter, *dataset);
        double scoreVal = scorer->test(*fitter, *validationData);

        if (verbose) {
          std::ostringstream out;
          out << "Score on batch: " << scoreTrain << std::endl
              << "Score on validation data: " << scoreVal;
          print(out);
        }

        visualizer->runVisualization(*fitter, *dataSource, fold, iteration);
        // Refine the model if neccessary
        monitor->pushToBuffer(numInstances, scoreVal, scoreTrain);
        size_t refinements = monitor->refinementsNecessary();
        while (refinements--) {
          fitter->adapt();
        }

        if (verbose) {
          std::ostringstream out;
          out << "###############"
              << "Iteration finished.";
          print(out);
        }
        iteration++;
      }
    }
    // Evaluate the final score on the validation data
    dataSource->reset();
    Dataset* validationData = dataSource->getValidationData();
    scores.push_back(scorer->test(*fitter, *validationData));
    delete monitor;  // release memory
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

  std::ostringstream out;
  out << "###############" << std::endl
      << "Mean score: " << meanScore << std::endl
      << "Standard deviation: " << stdDeviation;
  print(out);
  return meanScore;
}
} /* namespace datadriven */
} /* namespace sgpp */
