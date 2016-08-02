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
  std::vector<size_t> randomizedIndices(dataset.getNumberInstances());
  for (size_t i = 0; i < dataset.getNumberInstances(); i++) {
    randomizedIndices[i] = i;
  }

  //  for (size_t i : randomizedIndices) {
  //    std::cout << i << ", ";
  //  }
  std::cout << std::endl;

  shuffling->shuffle(randomizedIndices);

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
    size_t test_begin = fold * (dataset.getNumberInstances() / foldNumber);
    size_t test_end = test_begin + testSize;

    std::cout << "starting fold " << fold << " with full set:" << dataset.getNumberInstances()
              << ", test size:" << testSize << ", train size:" << trainSize << std::endl;
    std::cout << "test begins at:" << test_begin << ", test ends at:" << test_end << std::endl;
    // create testing & training datasets;
    auto testDataset = std::make_unique<Dataset>(testSize, dim);
    auto trainDataset = std::make_unique<Dataset>(trainSize, dim);

    // fill them
    DataVector tmpRow(dim);
    double tmpEntry = 0;
    size_t targetIdx = 0;
    // before test portion

    // std::cout << std::endl << "filling training portion" << std::endl << std::endl;

    for (size_t i = 0; i < test_begin; i++) {
      //      std::cout << "i=" << i << " with random index:" << randomizedIndices[i] << std::endl;
      dataset.getData().getRow(randomizedIndices[i], tmpRow);
      trainDataset->getData().setRow(i, tmpRow);
      tmpEntry = dataset.getTargets().get(randomizedIndices[i]);
      trainDataset->getTargets().set(i, tmpEntry);
    }

    //    std::cout << std::endl << "filling test portion" << std::endl << std::endl;

    // test portion
    targetIdx = 0;
    for (size_t i = test_begin; i < test_end; i++) {
      // std::cout << "i=" << targetIdx << " with random index:" << randomizedIndices[i] <<
      // std::endl;
      dataset.getData().getRow(randomizedIndices[i], tmpRow);
      testDataset->getData().setRow(targetIdx, tmpRow);
      tmpEntry = dataset.getTargets().get(randomizedIndices[i]);
      testDataset->getTargets().set(targetIdx, tmpEntry);
      targetIdx++;
    }

    //    std::cout << std::endl << "continue training portion" << std::endl << std::endl;

    // after test portion
    targetIdx = test_begin;
    for (size_t i = test_end; i < dataset.getNumberInstances(); i++) {
      //      std::cout << "i=" << targetIdx << " with random index:" << randomizedIndices[i] <<
      //      std::endl;
      dataset.getData().getRow(randomizedIndices[i], tmpRow);
      trainDataset->getData().setRow(targetIdx, tmpRow);
      tmpEntry = dataset.getTargets().get(randomizedIndices[i]);
      trainDataset->getTargets().set(targetIdx, tmpEntry);
      targetIdx++;
    }

    // fit model
    std::cout << "fitting model" << std::endl;
    model.fit(*trainDataset);

    // calculate testing error
    auto predictedValues = model.evaluate(testDataset->getData());

    // set score
    scores[fold] = (*metric)(*predictedValues, testDataset->getTargets());
    std::cout << "accuracy of fit:" << scores[fold] << std::endl;
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
