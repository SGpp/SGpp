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
    : Scorer{metric, shuffling, seed}, trainPortion{trainPortion} {}

Scorer* SplittingScorer::clone() const { return new SplittingScorer{*this}; }

// TODO(lettrich) :recycle
double SplittingScorer::calculateScore(ModelFittingBase& model, Dataset& dataset,
                                       double* stdDeviation) {
  // perform randomization of indices
  std::vector<size_t> randomizedIndices(dataset.getNumberInstances());
  randomizeIndices(dataset, randomizedIndices);

  // calculate size of testing and training portions
  size_t trainSize = std::lround(static_cast<double>(dataset.getNumberInstances()) * trainPortion);
  size_t testSize = dataset.getNumberInstances() - trainSize;
  size_t dim = dataset.getDimension();

  std::cout << "starting training with testing.\n"
            << "test size:" << testSize << "\n"
            << "train size:" << trainSize << "\n";

  // create test and train datasets.
  Dataset testDataset{testSize, dim};
  Dataset trainDataset{trainSize, dim};
  splitSet(dataset, trainDataset, testDataset, randomizedIndices);

  // fit model
  double score = train(model, trainDataset, testDataset);

  // refine it
  score = refine(model, testDataset);

  if (stdDeviation) {
    *stdDeviation = 0;
  }

  return score;
}

} /* namespace datadriven */
} /* namespace sgpp */
