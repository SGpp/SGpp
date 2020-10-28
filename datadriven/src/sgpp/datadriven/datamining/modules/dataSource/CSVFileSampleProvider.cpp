// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/dataSource/CSVFileSampleProvider.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/exception/file_exception.hpp>
#include <sgpp/datadriven/tools/CSVTools.hpp>

#include <string>
#include <vector>

namespace sgpp {

namespace datadriven {

CSVFileSampleProvider::CSVFileSampleProvider(DataShufflingFunctor* shuffling)
    : shuffling{shuffling}, dataset(Dataset{}), counter(0) {}

SampleProvider* CSVFileSampleProvider::clone() const {
  return dynamic_cast<SampleProvider*>(new CSVFileSampleProvider{*this});
}

size_t CSVFileSampleProvider::getDim() const {
  if (dataset.getDimension() != 0) {
    return dataset.getDimension();
  } else {
    throw base::file_exception{"No dataset loaded."};
  }
}

size_t CSVFileSampleProvider::getNumSamples() const {
  if (dataset.getDimension() != 0) {
    return dataset.getNumberInstances();
  } else {
    throw base::file_exception{"No dataset loaded."};
  }
}

void CSVFileSampleProvider::readFile(const std::string& fileName, bool hasTargets,
                                     size_t readinCutoff, std::vector<size_t> readinColumns,
                                     std::vector<double> readinClasses) {
  try {
    // call readCSV with skipfirstline set to true
    dataset = CSVTools::readCSVFromFile(fileName, true, hasTargets, readinCutoff, readinColumns,
                                        readinClasses);
  } catch (...) {
    // TODO(lettrich): catching all exceptions is bad design. Replace call to CSVTools with
    // exception safe implementation.
    throw base::data_exception{"Failed to parse CSV File."};
  }
}

Dataset* CSVFileSampleProvider::getNextSamples(size_t howMany) {
  if (dataset.getDimension() != 0) {
    return splitDataset(howMany);
  } else {
    throw base::file_exception("No dataset loaded.");
  }
}

Dataset* CSVFileSampleProvider::getAllSamples() {
  if (dataset.getDimension() != 0) {
    return this->getNextSamples(dataset.getNumberInstances());
  } else {
    throw base::file_exception{"No dataset loaded."};
  }
}

void CSVFileSampleProvider::readString(const std::string& input, bool hasTargets,
                                       size_t readinCutoff, std::vector<size_t> readinColumns,
                                       std::vector<double> readinClasses) {
  // try {
  //   dataset = CSVTools::readCSVFromString(input);
  // } catch (...) {
  // TODO(lettrich): catching all exceptions is bad design. Replace call to CSVTools with
  // exception safe implementation.
  throw base::data_exception{"Failed to parse CSV data. Reading from string is not implemented."};
  // }
}

Dataset* CSVFileSampleProvider::splitDataset(size_t howMany) {
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

void CSVFileSampleProvider::reset() { counter = 0; }

} /* namespace datadriven */
} /* namespace sgpp */
