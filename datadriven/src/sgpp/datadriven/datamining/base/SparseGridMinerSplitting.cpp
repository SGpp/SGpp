// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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

  // optimize lambda
  if (fitter->getFitterConfiguration().getRegularizationConfig().optimizeLambda_) {
    if (verbose) {
      std::ostringstream out;
      out << "Optimizing lambda";
      print(out);
    }

    double lambda = this->optimizeLambda(verbose);
    this->fitter->getFitterConfiguration().getRegularizationConfig().lambda_ = lambda;
    this->fitter->updateRegularization(lambda);

    if (verbose) {
      std::ostringstream out;
      out << "Final lambda: " << lambda;
      print(out);
    }
  }

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
            << "Iteration #" << (iteration) << std::endl
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

      visualizer->runVisualization(*fitter, *dataSource, 0, iteration);

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
}

double SparseGridMinerSplitting::optimizeLambda(bool verbose) {
  // 1 / phi
  double phiInversed = (std::sqrt(5) - 1) / 2;

  // 1 / phi^2
  double phiInversedSquared = (3 - std::sqrt(5)) / 2;

  // get configuration parameters
  double tolerance = fitter->getFitterConfiguration().getRegularizationConfig().optimizerTolerance_;
  double convergenceThreshold =
      fitter->getFitterConfiguration().getRegularizationConfig().convergenceThreshold_;
  double a = fitter->getFitterConfiguration().getRegularizationConfig().intervalA_;
  double b = fitter->getFitterConfiguration().getRegularizationConfig().intervalB_;

  double interval = std::abs(a - b);

  double c = a + phiInversedSquared * interval;
  double d = a + phiInversed * interval;

  // setup lambda for the first run
  fitter->getFitterConfiguration().getRegularizationConfig().lambda_ = c;

  double valueC = evaluateLambda(c, verbose);
  double valueD = evaluateLambda(d, verbose);
  double prevValueC = valueC;
  double prevValueD = valueD;

  int optimizerIteration = 0;

  while (interval > tolerance) {
    if (valueC < valueD) {
      b = d;
      d = c;
      valueD = valueC;
      interval = interval * phiInversed;
      c = a + phiInversedSquared * interval;
      valueC = evaluateLambda(c, verbose);
    } else {
      a = c;
      c = d;
      valueC = valueD;
      interval = interval * phiInversed;
      d = a + phiInversed * interval;
      valueD = evaluateLambda(d, verbose);
    }

    if (verbose) {
      std::ostringstream out;
      out << "###############"
          << "Lambda optimizer iteration #" << (optimizerIteration) << std::endl
          << "a: " << a << " b: " << b << " c: " << c << " d: " << d;
      print(out);
    }
    optimizerIteration++;

    // check for convergence
    if (valueC < valueD) {
      if (std::abs(valueC - prevValueC) < convergenceThreshold) {
        return a;
      }
    } else {
      if (std::abs(valueD - prevValueD) < convergenceThreshold) {
        return c;
      }
    }
    prevValueC = valueC;
    prevValueD = valueD;
  }

  if (valueC < valueD) {
    return a;
  }
  return b;
}

double SparseGridMinerSplitting::evaluateLambda(double lambda, bool verbose) {
  // update lambda of model (this has no effect if this is the first run)
  fitter->updateRegularization(lambda);

  // train the model (update loop)
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

    // Train model on new batch
    fitter->update(*dataset);
    iteration++;
  }

  // score the model
  double scoreVal = scorer->test(*fitter, *(dataSource->getValidationData()), true);

  // reset the fitter (but only the online part)
  fitter->resetTraining();

  return scoreVal;
}
}  // namespace datadriven
}  // namespace sgpp
