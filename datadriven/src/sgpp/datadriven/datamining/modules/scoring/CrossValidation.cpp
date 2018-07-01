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

void CrossValidation::calculateScore(ModelFittingBase& model, Dataset& dataset, double *scoreTrain,
                                       double *scoreTest, double* stdDeviation) {
  std::vector<double> scoresTrain(foldNumber);
  std::vector<double> scoresTest(foldNumber);

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
    train(model, trainDataset, testDataset, scoreTrain, scoreTest);
    scoresTrain[fold] = *scoreTrain;
    scoresTest[fold] = *scoreTest;
  }

  // calculate final score as AVG
  *scoreTrain = 0.0;
  *scoreTest = 0.0;
  for (size_t i = 0; i < foldNumber; i++) {
    *scoreTrain += scoresTrain[i];
    *scoreTest += scoresTest[i];
  }
  *scoreTrain /= static_cast<double>(foldNumber);
  *scoreTest /= static_cast<double>(foldNumber);

  // calculate std deviation if desired
  if (stdDeviation) {
    for (auto i : scoresTest) {
      *stdDeviation += std::pow(i - *scoreTest, 2);
    }
    *stdDeviation = sqrt((*stdDeviation) / (static_cast<double>(foldNumber) - 1.0));
  }
}
} /* namespace datadriven */
} /* namespace sgpp */
