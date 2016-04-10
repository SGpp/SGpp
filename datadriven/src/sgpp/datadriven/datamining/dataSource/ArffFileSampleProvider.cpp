/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ARFFWrapper.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun, Michael Lettrich
 */

#include <sgpp/datadriven/datamining/dataSource/ArffFileSampleProvider.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/base/exception/file_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <memory>
#include <string>

namespace sgpp {

namespace datadriven {

ArffFileSampleProvider::ArffFileSampleProvider() : dataset(Dataset()), counter(0) {}

ArffFileSampleProvider::~ArffFileSampleProvider() {}

size_t ArffFileSampleProvider::getDim() {
  if (dataset.getDimension() != 0) {
    return dataset.getDimension();
  } else {
    throw base::file_exception("No dataset loaded.");
  }
}

size_t ArffFileSampleProvider::getDatasetSize() {
  if (dataset.getDimension() != 0) {
    return dataset.getNumberInstances();
  } else {
    throw base::file_exception("No dataset loaded.");
  }
}

void ArffFileSampleProvider::readFile(const std::string& fileName) {
  dataset = ARFFTools::readARFF(fileName);
}

std::unique_ptr<Dataset> ArffFileSampleProvider::getNextSamples(size_t howMany) {
  if ((dataset.getDimension() != 0) && (counter + howMany) < dataset.getNumberInstances()) {
    return splitDataset(howMany);
  } else {
    throw base::data_exception("Demanded more samples then available.");
  }
}

std::unique_ptr<Dataset> ArffFileSampleProvider::getAllSamples() {
  if (dataset.getDimension() != 0) {
    auto tmpDataset = std::make_unique<Dataset>();
    (*tmpDataset) = dataset;
    return tmpDataset;
  } else {
    throw base::file_exception("No dataset loaded.");
  }
}

void ArffFileSampleProvider::readString(const std::string& input) {
  dataset = ARFFTools::readARFFFromString(input);
}

std::unique_ptr<Dataset> ArffFileSampleProvider::splitDataset(size_t howMany) {
  auto tmpDataset = std::make_unique<Dataset>(howMany, dataset.getDimension());

  base::DataMatrix& srcSamples = dataset.getData();
  base::DataVector& srcTargets = dataset.getTargets();

  base::DataMatrix& destSamples = tmpDataset->getData();
  base::DataVector& destTargets = tmpDataset->getTargets();

  base::DataVector tmpRow(srcSamples.getNcols());
  for (size_t i = counter; i < counter + howMany; ++i) {
    srcSamples.getRow(i, tmpRow);
    destSamples.setRow(i - counter, tmpRow);

    destTargets[i - counter] = srcTargets[i];
  }
  counter = counter + howMany;

  return tmpDataset;
}

} /* namespace datadriven */
} /* namespace sgpp */
