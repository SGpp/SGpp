/*
 * CrossValidationScorer.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#include <iostream>
#include <random>
#include <ctime>

#include <sgpp/globaldef.hpp>
#include "SimpleSplittingScorer.hpp"

using namespace SGPP::base;

namespace SGPP {
namespace datadriven {

SimpleSplittingScorer::SimpleSplittingScorer(Dataset& dataset, Metric* metric,
    double trainPortion, int seed): Scorer(metric),
  trainPortion(trainPortion), seed(seed) {
	if(this->seed){
		this->seed = std::time(NULL);
	}

}

SimpleSplittingScorer::~SimpleSplittingScorer() {
  // TODO Auto-generated destructor stub
}

void SimpleSplittingScorer::splitset(const DataMatrix& dataset,
                                     const DataVector& datasetValues, double trainPortion,
                                     DataMatrix& trainingSet,
									 DataVector& trainingValues,
									 DataMatrix& testSet,
									 DataVector& testValues, bool permute) {
  size_t dim = dataset.getNcols();
  size_t datasetSize = dataset.getNrows();  // size of data
  size_t trainSize = static_cast<size_t>(static_cast<double>(datasetSize) * trainPortion);
  size_t testSize = datasetSize - trainSize;

  // generate range of all indices and permute them if required
  std::vector<size_t> indices(datasetSize);

  for (size_t i = 0; i < datasetSize; i++ ) {
    indices.push_back(i);
  }

  if (permute == true) {
    std::mt19937 rndGen(seed);
    std::shuffle(indices.begin(), indices.end(), rndGen);
  }

  // fill data
  // first initialize data structures
  DataVector tmpLine(dim);
  trainingSet.resize(trainSize, dim);
  trainingValues.resize(trainSize);
  testSet.resize(testSize, dim);
  testValues.resize(testSize);

  // fill training data
  for (size_t j = 0; j < trainSize; j++) {
    size_t idx = indices[j];
    dataset.getRow(idx, tmpLine);

    for (size_t d = 0; d < dim; d++) {
      trainingSet.set(j, d, tmpLine[d]);
    }

    trainingValues.set(j, tmpLine[dim]);
  }

  // fill test data
  for (size_t j = trainSize; j < datasetSize; j++) {
    size_t idx = indices[j];
    dataset.getRow(idx, tmpLine);

    for (size_t d = 0; d < dim; d++) {
      testSet.set(j - trainSize, d, tmpLine[d]);
    }

    testValues.set(j - trainSize, tmpLine[dim]);
  }

}


} /* namespace datadriven */
} /* namespace SGPP */
