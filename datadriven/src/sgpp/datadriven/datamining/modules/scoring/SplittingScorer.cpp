/* Copyright (C) 2008-today The SG++ project

 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SplittingScorer.cpp
 *
 *  Created on:	07.10.2016
 *      Author: Michael Lettrich
 */
#include <sgpp/datadriven/datamining/modules/scoring/SplittingScorer.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

using base::DataVector;

SplittingScorer::SplittingScorer(Metric* metric, ShufflingFunctor* shuffling, int64_t seed,
                                 double trainPortion)
    : Scorer(metric, shuffling, seed), trainPortion(trainPortion) {}

SplittingScorer::~SplittingScorer() {}

double SplittingScorer::calculateScore(ModelFittingBase& model, Dataset& dataset,
                                       double* stdDeviation) {
  // perform randomization of indices
  std::vector<size_t> randomizedIndices(dataset.getNumberInstances());
  randomizeIndices(randomizedIndices, dataset.getNumberInstances());

  // calculate size of testing and training portions
  size_t trainSize = std::lround(static_cast<double>(dataset.getNumberInstances()) * trainPortion);
  size_t testSize = dataset.getNumberInstances() - trainSize;
  size_t dim = dataset.getDimension();

  std::cout << "starting training with testing. " << std::endl
            << "test size:" << testSize << std::endl
            << "train size:" << trainSize << std::endl;

  // create test and train datasets.
  auto testDataset = std::make_unique<Dataset>(testSize, dim);
  auto trainDataset = std::make_unique<Dataset>(trainSize, dim);
  splitSet(dataset, *trainDataset, *testDataset, trainSize, testSize, randomizedIndices);

  // fit model
  std::cout << "###############" << std::endl << "fitting model" << std::endl;
  double score = train(model, *trainDataset, *testDataset);
  std::cout << "###############" << std::endl
            << "accuracy of fit:" << score << std::endl
            << std::endl;

  if (stdDeviation) {
    *stdDeviation = 0;
  }

  return score;
}

} /* namespace datadriven */
} /* namespace sgpp */
