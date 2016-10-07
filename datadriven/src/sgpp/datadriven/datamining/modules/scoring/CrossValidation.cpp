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
    : metric(std::shared_ptr<Metric>(metric)),
      shuffling(std::shared_ptr<ShufflingFunctor>(shuffling)),
      foldNumber(foldNumber) {
  if (seed != -1) {
    shuffling->setSeed(seed);
  }
}

double CrossValidation::calculateScore(ModelFittingBase& model, Dataset& dataset,
                                       double* stdDeviation) {
  std::vector<double> scores(foldNumber);

  // perform randomization of indices
  std::vector<size_t> randomizedIndices(dataset.getNumberInstances());
  for (size_t i = 0; i < dataset.getNumberInstances(); i++) {
    randomizedIndices[i] = i;
  }

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

    for (size_t i = 0; i < test_begin; i++) {
      dataset.getData().getRow(randomizedIndices[i], tmpRow);
      trainDataset->getData().setRow(i, tmpRow);
      tmpEntry = dataset.getTargets().get(randomizedIndices[i]);
      trainDataset->getTargets().set(i, tmpEntry);
    }

    // test portion
    targetIdx = 0;
    for (size_t i = test_begin; i < test_end; i++) {
      dataset.getData().getRow(randomizedIndices[i], tmpRow);
      testDataset->getData().setRow(targetIdx, tmpRow);
      tmpEntry = dataset.getTargets().get(randomizedIndices[i]);
      testDataset->getTargets().set(targetIdx, tmpEntry);
      targetIdx++;
    }

    // after test portion
    targetIdx = test_begin;
    for (size_t i = test_end; i < dataset.getNumberInstances(); i++) {
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
    // auto predictedValues = model.evaluate(testDataset->getData());

    auto predictedValues = std::make_unique<DataVector>(testDataset->getNumberInstances());
    model.evaluate(testDataset->getData(), *predictedValues);
    // set score
    scores[fold] = (*metric)(*predictedValues, testDataset->getTargets());
    std::cout << "accuracy of fit:" << scores[fold] << std::endl;
    predictedValues.reset();
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
