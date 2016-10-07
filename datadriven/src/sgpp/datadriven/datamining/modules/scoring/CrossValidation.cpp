/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CrossValidation.cpp
 *
 *  Created on: 31.07.2016
 *      Author: Michael Lettrich
 */

#include "CrossValidation.hpp"
#include <iostream>
#include <numeric>
#include <vector>

namespace sgpp {
namespace datadriven {

CrossValidation::~CrossValidation() {}

CrossValidation::CrossValidation(Metric* metric, ShufflingFunctor* shuffling, int64_t seed,
                                 size_t foldNumber)
    : Scorer(metric, shuffling), foldNumber(foldNumber) {}

double CrossValidation::calculateScore(ModelFittingBase& model, Dataset& dataset,
                                       double* stdDeviation) {
  std::vector<double> scores(foldNumber);

  // perform randomization of indices
  std::vector<size_t> randomizedIndices(dataset.getNumberInstances());
  randomizeIndices(randomizedIndices, dataset.getNumberInstances());

  // perform actual folding

  for (size_t fold = 0; fold < foldNumber; fold++) {
    // calculate size of the fold
    size_t testSize = dataset.getNumberInstances() / foldNumber;
    // for last fold add all remaining items
    if (fold == foldNumber - 1) {
      testSize += dataset.getNumberInstances() % foldNumber;
    }
    size_t trainSize = dataset.getNumberInstances() - testSize;
    size_t dim = dataset.getDimension();

    // calculate where the training portion starts
    size_t offset = fold * (dataset.getNumberInstances() / foldNumber);

    std::cout << "starting fold " << fold << " with full set:" << dataset.getNumberInstances()
              << ", test size:" << testSize << ", train size:" << trainSize << std::endl;
    // create testing & training datasets;
    auto testDataset = std::make_unique<Dataset>(testSize, dim);
    auto trainDataset = std::make_unique<Dataset>(trainSize, dim);

    // fill them
    splitSet(dataset, *trainDataset, *testDataset, trainSize, testSize, randomizedIndices, offset);

    // fit model
    std::cout << "fitting model" << std::endl;
    scores[fold] = train(model, *trainDataset, *testDataset);
    std::cout << "accuracy of fit:" << scores[fold] << std::endl;
  }

  // calculate final score as AVG
  double avgScore = 0;
  for (auto i : scores) {
    avgScore += i;
  }
  avgScore = avgScore / static_cast<double>(foldNumber);

  // calculate std deviation if desired
  if (stdDeviation) {
    for (auto i : scores) {
      *stdDeviation += std::pow(i - avgScore, 2);
    }
    *stdDeviation = sqrt((*stdDeviation) / (static_cast<double>(foldNumber) - 1.0));
  }

  return avgScore;
}
} /* namespace datadriven */
} /* namespace sgpp */
