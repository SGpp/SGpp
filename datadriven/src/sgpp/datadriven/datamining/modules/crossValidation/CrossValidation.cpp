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
#include <numeric>
#include <vector>

namespace sgpp {
namespace datadriven {

CrossValidation::~CrossValidation() {}

CrossValidation::CrossValidation(std::shared_ptr<Metric> metric,
                                 std::shared_ptr<ShufflingFunctor> shuffling, int64_t seed)
    : metric(metric), shuffling(shuffling) {
  if (seed != -1) {
    shuffling->setSeed(seed);
  }
}

double CrossValidation::calculateScore(ModelFittingBase& model, Dataset& dataset, size_t foldNumber,
                                       std::shared_ptr<double> stdDeviation) {
  std::vector<double> scores(foldNumber);

  // perform randomization of indices
  std::vector<size_t> randomizedIndices;
  for (size_t i = 0; i < dataset.getNumberInstances(); i++) {
    randomizedIndices.emplace_back(i);
  }
  shuffling->shuffle(randomizedIndices);

  // perform actual folding
  for (size_t fold = 0; fold < foldNumber; fold++) {
    // calculate size of the fold
    size_t testSize = dataset.getNumberInstances() / foldNumber;
    if (fold == foldNumber - 1) {
      testSize += dataset.getNumberInstances() % foldNumber;
    }
    size_t trainSize = dataset.getNumberInstances() - testSize;
    size_t dim = dataset.getDimension();

    // calculate where the training portion starts
    size_t test_begin = fold * (dataset.getNumberInstances() / foldNumber);
    size_t test_end = test_begin + testSize;

    // create testing & training datasets;
    auto testDataset = std::make_unique<Dataset>(testSize, dim);
    auto trainDataset = std::make_unique<Dataset>(trainSize, dim);

    // fill them
    DataVector tmpRow(dim);
    double tmpEntry = 0;

    // before test portion
    for (size_t i = 0; i < test_begin; i++) {
      dataset.getData().getRow(randomizedIndices[i], tmpRow);
      trainDataset->getData().setRow(i, tmpRow);
      tmpEntry = dataset.getTargets().get(randomizedIndices[i]);
      trainDataset->getTargets().set(i, tmpEntry);
    }

    // test portion
    for (size_t i = test_begin; i < test_end; i++) {
      dataset.getData().getRow(randomizedIndices[i], tmpRow);
      testDataset->getData().setRow(i, tmpRow);
      tmpEntry = dataset.getTargets().get(randomizedIndices[i]);
      testDataset->getTargets().set(i, tmpEntry);
    }

    // after test portion
    for (size_t i = test_end; i < dataset.getNumberInstances(); i++) {
      dataset.getData().getRow(randomizedIndices[i], tmpRow);
      trainDataset->getData().setRow(i, tmpRow);
      tmpEntry = dataset.getTargets().get(randomizedIndices[i]);
      trainDataset->getTargets().set(i, tmpEntry);
    }

    // fit model
    model.fit(*trainDataset);

    // calculate testing error
    auto predictedValues = model.evaluate(testDataset->getData());

    // set score
    scores.emplace_back((*metric)(*predictedValues, testDataset->getTargets()));
    predictedValues.release();
  }

  // calculate final score as AVG
  double avgScore = std::accumulate(scores.begin(), scores.end(), 0);
  avgScore /= static_cast<double>(foldNumber);

  // calculate std deviation if desired
  if (stdDeviation) {
    *stdDeviation = std::accumulate(scores.begin(), scores.end(), 0,
                                    [avgScore](double a, double b) { return a + (b - avgScore); });
    *stdDeviation = sqrt((*stdDeviation) / (static_cast<double>(foldNumber) - 1.0));
  }

  return avgScore;
}
} /* namespace datadriven */
} /* namespace sgpp */
