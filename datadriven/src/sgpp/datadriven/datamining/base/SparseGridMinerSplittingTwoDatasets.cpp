// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/base/SparseGridMinerSplittingTwoDatasets.hpp>

#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorFactory.hpp>

#include <iostream>
#include <vector>

namespace sgpp {
namespace datadriven {

SparseGridMinerSplittingTwoDatasets::SparseGridMinerSplittingTwoDatasets(
    std::vector<DataSourceSplitting*> dataSource, ModelFittingBase* fitter, Scorer* scorer,
    Visualizer* visualizer)
    : SparseGridMiner(fitter, scorer, visualizer),
      dataSourceP{dataSource[0]},
      dataSourceQ{dataSource[1]} {}

double SparseGridMinerSplittingTwoDatasets::learn(bool verbose) {
  fitter->verboseSolver = verbose;
  // Setup refinement monitor
  RefinementMonitorFactory monitorFactory;
  RefinementMonitor* monitor = monitorFactory.createRefinementMonitor(
      fitter->getFitterConfiguration().getRefinementConfig());

  // We use only the parameters set for the first dataSource instance; Batching is not implemented
  for (size_t epoch = 0; epoch < dataSourceP->getConfig().epochs_; epoch++) {
    if (verbose) {
      std::cout << "###############"
                << "Starting training epoch #" << epoch << std::endl;
    }
    dataSourceP->reset();
    dataSourceQ->reset();
    // Process dataset iteratively
    size_t iteration = 0;
    while (true) {
      std::unique_ptr<Dataset> datasetP(dataSourceP->getNextSamples());
      std::unique_ptr<Dataset> datasetQ(dataSourceQ->getNextSamples());

      size_t numInstancesP = datasetP->getNumberInstances();
      size_t numInstancesQ = datasetQ->getNumberInstances();

      if (numInstancesP == 0 || numInstancesQ == 0) {
        // One of the sources does not provide any more samples
        break;
      }
      if (verbose) {
        std::cout << "###############"
                  << "Itertation #" << (iteration++) << std::endl
                  << "Batch size: " << numInstancesP << ", " << numInstancesQ << std::endl;
      }
      // Train model on new batch
      fitter->update(*datasetP, *datasetQ);

      // Evaluate the score on the training and validation data, for each of the input datasets
      // We average the scores of the two datasets
      double avgScoreTrain =
          (scorer->test(*fitter, *datasetP) + scorer->test(*fitter, *datasetQ)) / 2.;
      double avgScoreVal = (scorer->test(*fitter, *(dataSourceP->getValidationData())) +
                            scorer->test(*fitter, *(dataSourceQ->getValidationData()))) /
                           2.;

      if (verbose) {
        std::cout << "Score on batch: " << avgScoreTrain << std::endl
                  << "Score on validation data: " << avgScoreVal << std::endl;
      }
      // Refine the model if neccessary
      monitor->pushToBuffer(numInstancesP + numInstancesQ, avgScoreVal, avgScoreTrain);
      size_t refinements = monitor->refinementsNecessary();
      while (refinements--) {
        fitter->adapt();
      }
      if (verbose) {
        std::cout << "###############"
                  << "Iteration finished." << std::endl;
      }
    }
  }
  delete monitor;  // release memory
  return (scorer->test(*fitter, *(dataSourceP->getValidationData())) +
          scorer->test(*fitter, *(dataSourceQ->getValidationData()))) /
         2.;
}
} /* namespace datadriven */
} /* namespace sgpp */
