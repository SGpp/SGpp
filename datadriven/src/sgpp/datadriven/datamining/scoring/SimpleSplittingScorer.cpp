/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CrossValidationScorer.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#include <random>
#include <ctime>
#include <iostream>

#include <vector>

#include <sgpp/globaldef.hpp>                                    //NOLINT (is no systemheader)
#include <sgpp/datadriven/datamining/SimpleSplittingScorer.hpp>  //NOLINT (is no systemheader)
#include <sgpp/base/tools/json/json_exception.hpp>               //NOLINT (is no systemheader)

using namespace sgpp::base;  // NOLINT

namespace sgpp {
namespace datadriven {

SimpleSplittingScorer::SimpleSplittingScorer(std::shared_ptr<Metric> metric,
                                             std::shared_ptr<ModelFittingBase> fitter,
                                             datadriven::DataMiningConfiguration config)
    : Scorer(metric, fitter) {
  try {
    trainPortion = config["trainPortion"].getDouble();
    // TODO(lettrich): if seed not set use random seed
    seed = config["seed"].getUInt();
  } catch (json::json_exception& e) {
    std::cout << e.what() << std::endl;
  }
}

SimpleSplittingScorer::~SimpleSplittingScorer() {
  // TODO(lettrich): Auto-generated destructor stub
}

double SimpleSplittingScorer::getScore(Dataset& dataset) {
  // initialize data structures
  size_t dim = dataset.getDimension();
  size_t datasetSize = dataset.getNumberInstances();  // size of data
  size_t trainSize = static_cast<size_t>(static_cast<double>(datasetSize) * trainPortion);
  size_t testSize = datasetSize - trainSize;
  Dataset trainingSet(trainSize, dim);
  Dataset testSet(testSize, dim);
  splitset(dataset, trainingSet, testSet);

  // fit the model to the training dataset
  fitter->update(trainingSet);

  // evaluate the model on the test set and return the metric
  DataVector predictedValues(testSize);
  fitter->evaluate(testSet.getData(), predictedValues);
  return (*metric)(predictedValues, testSet.getTargets());
}

void SimpleSplittingScorer::splitset(Dataset& dataset, Dataset& trainingSet, Dataset& testSet,
                                     bool permute) {
  size_t dim = dataset.getDimension();
  size_t datasetSize = dataset.getNumberInstances();  // size of data
  size_t trainSize = trainingSet.getNumberInstances();
  size_t testSize = testSet.getNumberInstances();

  // generate range of all indices and permute them if required
  std::vector<size_t> indices(datasetSize);

  for (size_t i = 0; i < datasetSize; i++) {
    indices.push_back(i);
  }

  if (permute == true) {
    std::mt19937 rndGen(seed);
    std::shuffle(indices.begin(), indices.end(), rndGen);
  }

  // fill data
  // first initialize data structures
  DataVector tmpLine(dim);
  trainingSet.getData().resize(trainSize, dim);
  trainingSet.getTargets().resize(trainSize);
  testSet.getData().resize(testSize, dim);
  testSet.getTargets().resize(testSize);

  // fill training data
  for (size_t j = 0; j < trainSize; j++) {
    size_t idx = indices[j];
    dataset.getData().getRow(idx, tmpLine);

    for (size_t d = 0; d < dim; d++) {
      trainingSet.getData().set(j, d, tmpLine[d]);
    }

    trainingSet.getTargets().set(j, tmpLine[dim]);
  }

  // fill test data
  for (size_t j = trainSize; j < datasetSize; j++) {
    size_t idx = indices[j];
    dataset.getData().getRow(idx, tmpLine);

    for (size_t d = 0; d < dim; d++) {
      testSet.getData().set(j - trainSize, d, tmpLine[d]);
    }

    testSet.getTargets().set(j - trainSize, tmpLine[dim]);
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
