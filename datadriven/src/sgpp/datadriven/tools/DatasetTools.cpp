// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/tools/DatasetTools.hpp>

#include <iostream>
#include <vector>

namespace sgpp {
namespace datadriven {

void DatasetTools::splitset(base::DataMatrix& dataset,
                            base::DataVector& datasetValues, size_t kFold,
                            std::vector<base::DataMatrix>& trainingSets,
                            std::vector<base::DataVector>& trainingSetValues,
                            std::vector<base::DataMatrix>& testSets,
                            std::vector<base::DataVector>& testSetValues, bool verbose) {
  trainingSets.clear();
  trainingSetValues.clear();
  testSets.clear();
  testSetValues.clear();

  size_t dim = dataset.getNcols();

  base::DataVector p(dim);
  base::DataVector tmp(dim);

  std::vector<size_t> partitionSizes(kFold);  // size of partition
  std::vector<size_t> partitionStartIndices(kFold + 1);  // index of partition
  size_t datasetSize = dataset.getNrows();  // size of data

  // set size of partitions
  if (verbose) {
    std::cout << "kfold partitions: ";
  }

  partitionStartIndices[0] = 0;

  for (size_t i = 0; i < kFold - 1; i++) {
    partitionSizes[i] = datasetSize / kFold;
    partitionStartIndices[i + 1] = partitionStartIndices[i] + partitionSizes[i];

    if (verbose) {
      std::cout << partitionSizes[i] << " ";
    }
  }

  partitionStartIndices[kFold] = datasetSize;
  partitionSizes[kFold - 1] = datasetSize - (kFold - 1) * (datasetSize / kFold);

  if (verbose) {
    std::cout << partitionSizes[kFold - 1] << std::endl;
  }

  if (verbose) {
    std::cout << "kfold indices: ";

    for (size_t i = 0; i <= kFold; i++) {
      std::cout << partitionStartIndices[i] << " ";
    }

    std::cout << std::endl;
  }

  // fill data
  for (size_t i = 0; i < kFold; i++) {
    trainingSets.emplace_back(dataset.getNrows() - partitionSizes[i],
                              dataset.getNcols());
    trainingSetValues.emplace_back(dataset.getNrows() - partitionSizes[i]);
    testSets.emplace_back(partitionSizes[i], dataset.getNcols());
    testSetValues.emplace_back(partitionSizes[i]);

    size_t local_test = 0;
    size_t local_train = 0;

    for (size_t j = 0; j < dataset.getNrows(); j++) {
      dataset.getRow(j, p);

      if (partitionStartIndices[i] <= j && j < partitionStartIndices[i + 1]) {
        testSets[i].setRow(local_test, p);
        testSetValues[i].set(local_test, datasetValues[j]);
        local_test++;
      } else {
        trainingSets[i].setRow(local_train, p);
        trainingSetValues[i].set(local_train, datasetValues[j]);
        local_train++;
      }
    }
  }
}

}  // namespace datadriven
}  // namespace sgpp

