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

#include <sgpp/datadriven/datamining/modules/scoring/CrossValidation.hpp>

#include <iostream>
#include <vector>

namespace sgpp {
namespace datadriven {

CrossValidation::CrossValidation(Metric* metric, ShufflingFunctor* shuffling, int64_t seed,
                                 size_t foldNumber)
    : Scorer{metric, shuffling}, foldNumber{foldNumber} {}

Scorer* CrossValidation::clone() const { return new CrossValidation{*this}; }

double CrossValidation::calculateScore(ModelFittingBase& model, Dataset& dataset,
                                       double* stdDeviation) {
  std::vector<double> scores(foldNumber);

  // perform randomization of indices
  std::vector<size_t> randomizedIndices(dataset.getNumberInstances());
  randomizeIndices(dataset, randomizedIndices);

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
              << ", test size:" << testSize << ", train size:" << trainSize << "\n";
    // create testing & training datasets;
    Dataset testDataset{testSize, dim};
    Dataset trainDataset{trainSize, dim};

    // fill them
    splitSet(dataset, trainDataset, testDataset, randomizedIndices, offset);

    // fit model
    scores[fold] = train(model, trainDataset, testDataset);

    scores[fold] = refine(model, testDataset);
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
