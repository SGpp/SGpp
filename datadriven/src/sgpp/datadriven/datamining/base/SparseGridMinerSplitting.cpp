/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseGridMinerSplitting.cpp
 *
 *  Created on: Jul 23, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/datamining/base/SparseGridMinerSplitting.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorFactory.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {


SparseGridMinerSplitting::SparseGridMinerSplitting(DataSourceSplitting* dataSource,
                                                   ModelFittingBase* fitter, Scorer* scorer,
                                                   Visualizer* visualizer)
    : SparseGridMiner(fitter, scorer, visualizer), dataSource{dataSource} {}

double SparseGridMinerSplitting::learn(bool verbose) {
#ifdef USE_SCALAPACK
  if (fitter->getFitterConfiguration().getParallelConfig().scalapackEnabled_) {
    auto processGrid = fitter->getProcessGrid();
    if (!processGrid->isProcessInGrid()) {
      return 0.0;
    }
  }
#endif /* USE_SCALAPACK */

  fitter->verboseSolver = verbose;
  // Setup refinement monitor
  RefinementMonitorFactory monitorFactory;
  RefinementMonitor* monitor = monitorFactory.createRefinementMonitor(
      fitter->getFitterConfiguration().getRefinementConfig());

  for (size_t epoch = 0; epoch < dataSource->getConfig().epochs; epoch++) {
    if (verbose) {
      std::ostringstream out;
      out << "###############"
          << "Starting training epoch #" << epoch;
      print(out);
    }
    dataSource->reset();
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
            << "Itertation #" << (iteration) << std::endl
            << "Batch size: " << numInstances;
        print(out);
      }
      // Train model on new batch
      fitter->update(*dataset);


      // Evaluate the score on the training and validation data
      double scoreTrain = scorer->test(*fitter, *dataset);
      double scoreVal = scorer->test(*fitter, *(dataSource->getValidationData()));


      if (verbose) {
        std::ostringstream out;
        out << "Score on batch: " << scoreTrain << std::endl
            << "Score on validation data: " << scoreVal;
        print(out);
      }

      visualizer->visualize(*fitter, 1, iteration);

      // Refine the model if neccessary
      monitor->pushToBuffer(numInstances, scoreVal, scoreTrain);
      size_t refinements = monitor->refinementsNecessary();
      while (refinements--) {
        fitter->refine();
      }
      if (verbose) {
        print("###############Iteration finished.");
      }
      iteration++;
    }
  }
  return scorer->test(*fitter, *(dataSource->getValidationData()));
}  // namespace datadriven
}  // namespace datadriven
}  // namespace sgpp
