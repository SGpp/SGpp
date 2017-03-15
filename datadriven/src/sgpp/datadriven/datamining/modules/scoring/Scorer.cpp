/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Scorer.cpp
 *
 * Created on: Oct 7, 2016
 *     Author: Michael Lettrich
 */

#include "Scorer.hpp"
#include <vector>

namespace sgpp {
namespace datadriven {

Scorer::Scorer(Metric* metric, ShufflingFunctor* shuffling, int64_t seed)
    : metric{std::unique_ptr<Metric>{metric}},
      shuffling{std::unique_ptr<ShufflingFunctor>{shuffling}} {
  if (seed != -1) {
    shuffling->setSeed(seed);
  }
}

Scorer::Scorer(const Scorer& rhs) {
  metric = std::unique_ptr<Metric>{rhs.metric->clone()};
  shuffling = std::unique_ptr<ShufflingFunctor>{rhs.shuffling->clone()};
}

Scorer& Scorer::operator=(const Scorer& rhs) {
  if (this != &rhs) {
    metric = std::unique_ptr<Metric>{rhs.metric->clone()};
    shuffling = std::unique_ptr<ShufflingFunctor>{rhs.shuffling->clone()};
  }
  return *this;
}

void Scorer::randomizeIndices(std::vector<size_t>& randomizedIndices) {
  for (size_t i = 0; i < randomizedIndices.size(); i++) {
    randomizedIndices[i] = i;
  }
  shuffling->shuffle(randomizedIndices);
}

void Scorer::splitSet(const Dataset& dataset, Dataset& trainDataset, Dataset& testDataset,
                      const std::vector<size_t>& randomizedIndices, size_t offset) {
  size_t testSize = testDataset.getNumberInstances();
  size_t dim = dataset.getDimension();
  DataVector tmpRow(dim);
  double tmpEntry = 0;
  size_t targetIdx = 0;

  size_t testBegin = offset;
  size_t testEnd = testBegin + testSize;

  // before test portion
  for (size_t i = 0; i < offset; i++) {
    dataset.getData().getRow(randomizedIndices[i], tmpRow);
    trainDataset.getData().setRow(i, tmpRow);
    tmpEntry = dataset.getTargets().get(randomizedIndices[i]);
    trainDataset.getTargets().set(i, tmpEntry);
  }

  // test portion
  targetIdx = 0;
  for (size_t i = offset; i < testEnd; i++) {
    dataset.getData().getRow(randomizedIndices[i], tmpRow);
    testDataset.getData().setRow(targetIdx, tmpRow);
    tmpEntry = dataset.getTargets().get(randomizedIndices[i]);
    testDataset.getTargets().set(targetIdx, tmpEntry);
    targetIdx++;
  }

  // after test portion
  targetIdx = offset;
  for (size_t i = testEnd; i < dataset.getNumberInstances(); i++) {
    dataset.getData().getRow(randomizedIndices[i], tmpRow);
    trainDataset.getData().setRow(targetIdx, tmpRow);
    tmpEntry = dataset.getTargets().get(randomizedIndices[i]);
    trainDataset.getTargets().set(targetIdx, tmpEntry);
    targetIdx++;
  }
}

double Scorer::train(ModelFittingBase& model, Dataset& trainDataset, Dataset& testDataset) {
  model.fit(trainDataset);

  DataVector predictedValues{testDataset.getNumberInstances()};
  model.evaluate(testDataset.getData(), predictedValues);
  // set score
  return metric->measure(predictedValues, testDataset.getTargets());
}

} /* namespace datadriven */
} /* namespace sgpp */
