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

// FOR EVALUATION ONLY
#include <sgpp/datadriven/datamining/base/evaluationtools.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorFactory.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerSplitting.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <iostream>

using evalu::getTime;

namespace sgpp {
namespace datadriven {

SparseGridMinerSplitting::SparseGridMinerSplitting(DataSourceSplitting* dataSource,
                                                   ModelFittingBase* fitter, Scorer* scorer)
    : SparseGridMiner(fitter, scorer), dataSource{dataSource} {}

double SparseGridMinerSplitting::learn(bool verbose) {
  fitter->verboseSolver = verbose;
  // Setup refinement monitor
  RefinementMonitorFactory monitorFactory;
  RefinementMonitor* monitor = monitorFactory.createRefinementMonitor(
      fitter->getFitterConfiguration().getRefinementConfig());
  for (size_t epoch = 0; epoch < dataSource->getConfig().epochs; epoch++) {
    if (verbose) {
      std::cout << "###############"
                << "Starting training epoch #" << epoch << std::endl;
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
        std::cout << "###############" << getTime() << "Itertation #" << (iteration++) << std::endl
                  << "Batch size: " << numInstances << std::endl;
      }
      // Train model on new batch
      std::cout << "TIME BEVOR UPDATE" << evalu::getMSek();
      fitter->update(*dataset);
      std::cout << "TIME AFTER UPDATE" << evalu::getMSek();

      // Evaluate the score on the training and validation data
      double scoreTrain = scorer->test(*fitter, *dataset);
      double scoreVal = scorer->test(*fitter, *(dataSource->getValidationData()));

      if (verbose) {
        std::cout << getTime() << std::endl
                  << "Score on batch: " << scoreTrain << std::endl
                  << "Score on validation data: " << scoreVal << std::endl;
      }
      // Refine the model if neccessary
      monitor->pushToBuffer(numInstances, scoreVal, scoreTrain);
      size_t refinements = monitor->refinementsNecessary();
      while (refinements--) {
        std::cout << getTime();
        fitter->refine();
        std::cout << getTime();
      }
      std::cout << std::endl << evalu::getMSek();
      if (verbose) {
        std::cout << "###############"
                  << "Iteration finished." << std::endl;
      }
    }
  }
  return scorer->test(*fitter, *(dataSource->getValidationData()));
}
} /* namespace datadriven */
} /* namespace sgpp */
