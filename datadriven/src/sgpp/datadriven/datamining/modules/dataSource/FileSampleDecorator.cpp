/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FileSampleDecorator.cpp
 *
 *  Created on: 12.12.2016
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleDecorator.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

FileSampleDecorator::FileSampleDecorator(FileSampleProvider *const fileSampleProvider)
    : FileSampleProvider{},
      fileSampleProvider(std::unique_ptr<FileSampleProvider>{fileSampleProvider}) {}

FileSampleDecorator::FileSampleDecorator(const FileSampleDecorator &rhs) : FileSampleProvider{} {
  fileSampleProvider = std::unique_ptr<FileSampleProvider>{
      dynamic_cast<FileSampleProvider *>(rhs.fileSampleProvider->clone())};
}

FileSampleDecorator &FileSampleDecorator::operator=(const FileSampleDecorator &rhs) {
  if (&rhs == this) {
    return *this;
  }
  fileSampleProvider = std::unique_ptr<FileSampleProvider>{
      dynamic_cast<FileSampleProvider *>(rhs.fileSampleProvider->clone())};
  return *this;
}

Dataset *FileSampleDecorator::getNextSamples(size_t howMany) {
  return fileSampleProvider->getNextSamples(howMany);
}

Dataset *FileSampleDecorator::getAllSamples() { return fileSampleProvider->getAllSamples(); }

size_t FileSampleDecorator::getDim() const { return fileSampleProvider->getDim(); }

size_t FileSampleDecorator::getNumSamples() const { return fileSampleProvider->getNumSamples(); }

void FileSampleDecorator::readFile(const std::string &fileName,
                                   bool hasTargets,
                                   size_t readinCutoff,
                                   std::vector<size_t> readinColumns,
                                   std::vector<double> readinClasses) {
  fileSampleProvider->readFile(fileName, hasTargets, readinCutoff, readinColumns, readinClasses);
}

void FileSampleDecorator::readString(const std::string &input,
                                     bool hasTargets,
                                     size_t readinCutoff,
                                     std::vector<size_t> readinColumns,
                                     std::vector<double> readinClasses) {
  fileSampleProvider->readString(input, hasTargets, readinCutoff, readinColumns, readinClasses);
}
} /* namespace datadriven */
} /* namespace sgpp */
