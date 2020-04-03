// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/exception/file_exception.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <string>
#include <vector>

namespace sgpp {

namespace datadriven {

ArffFileSampleProvider::ArffFileSampleProvider(DataShufflingFunctor* shuffling)
    : shuffling{shuffling}, dataset(Dataset{}), counter(0) {}

SampleProvider* ArffFileSampleProvider::clone() const {
  return dynamic_cast<SampleProvider*>(new ArffFileSampleProvider{*this});
}

size_t ArffFileSampleProvider::getDim() const {
  if (dataset.getDimension() != 0) {
    return dataset.getDimension();
  } else {
    throw base::file_exception{"No dataset loaded."};
  }
}

size_t ArffFileSampleProvider::getNumSamples() const {
  if (dataset.getDimension() != 0) {
    return dataset.getNumberInstances();
  } else {
    throw base::file_exception{"No dataset loaded."};
  }
}

void ArffFileSampleProvider::readFile(const std::string& fileName, bool hasTargets,
                                      size_t readinCutoff, std::vector<size_t> readinColumns,
                                      std::vector<double> readinClasses) {
  try {
    dataset = ARFFTools::readARFFFromFile(fileName, hasTargets, readinCutoff, readinColumns,
                                          readinClasses);
  } catch (...) {
    // TODO(lettrich): catching all exceptions is bad design. Replace call to ARFFTools with
    // exception safe implementation.
    throw base::data_exception{"Failed to parse ARFF File."};
  }
}

Dataset* ArffFileSampleProvider::getNextSamples(size_t howMany) {
  if (dataset.getDimension() != 0) {
    return splitDataset(howMany);
  } else {
    throw base::file_exception("No dataset loaded.");
  }
}

Dataset* ArffFileSampleProvider::getAllSamples() {
  if (dataset.getDimension() != 0) {
    return this->getNextSamples(dataset.getNumberInstances());
  } else {
    throw base::file_exception{"No dataset loaded."};
  }
}

void ArffFileSampleProvider::readString(const std::string& input, bool hasTargets,
                                        size_t readinCutoff, std::vector<size_t> readinColumns,
                                        std::vector<double> readinClasses) {
  try {
    dataset = ARFFTools::readARFFFromString(input, hasTargets, readinCutoff, readinColumns,
                                            readinClasses);
  } catch (...) {
    // TODO(lettrich): catching all exceptions is bad design. Replace call to ARFFTools with
    // exception safe implementation.
    throw base::data_exception{"Failed to parse ARFF data."};
  }
}

Dataset* ArffFileSampleProvider::splitDataset(size_t howMany) {
  const size_t size = counter + howMany <= dataset.getNumberInstances()
                          ? howMany
                          : dataset.getNumberInstances() - counter;
  auto tmpDataset = std::make_unique<Dataset>(size, dataset.getDimension());

  base::DataMatrix& srcSamples = dataset.getData();
  base::DataVector& srcTargets = dataset.getTargets();

  base::DataMatrix& destSamples = tmpDataset->getData();
  base::DataVector& destTargets = tmpDataset->getTargets();

  base::DataVector tmpRow(srcSamples.getNcols());

  // copy "size" rows beginning from "counter" to the new dataset.
  for (size_t i = counter; i < counter + size; ++i) {
    size_t srcIdx = shuffling != nullptr ? (*shuffling)(i, dataset.getNumberInstances()) : i;
    srcSamples.getRow(srcIdx, tmpRow);
    destSamples.setRow(i - counter, tmpRow);

    destTargets[i - counter] = srcTargets[srcIdx];
  }
  counter = counter + size;

  return tmpDataset.release();
}

void ArffFileSampleProvider::reset() { counter = 0; }

} /* namespace datadriven */
} /* namespace sgpp */
