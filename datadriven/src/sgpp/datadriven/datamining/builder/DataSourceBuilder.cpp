/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FileSampleProviderBuilder.cpp
 *
 *  Created on: 16.05.2016
 *      Author: Michael Lettrich
 */

#include "DataSourceBuilder.hpp"

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceState.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/GzipFileSampleDecorator.hpp>

#include <algorithm>
#include <string>

const std::string arff = "ARFF";

namespace sgpp {
namespace datadriven {

DataSourceBuilder::DataSourceBuilder()
    : fileType(NONE), filePath(""), isCompressed(false), batchSize(0), numBatches(0) {}

DataSourceBuilder::~DataSourceBuilder() {}

DataSourceBuilder& DataSourceBuilder::withFileType(std::string fileType) {
  if (arff.compare(fileType) == 0) {
    this->fileType = ARFF;
  } else {
    base::data_exception("Unknown file type");
  }
  return *this;
}

DataSourceBuilder& DataSourceBuilder::withCompression(bool isCompressed) {
  this->isCompressed = isCompressed;
  return *this;
}

DataSourceBuilder& DataSourceBuilder::withPath(std::string filePath) {
  this->filePath = filePath;
  return *this;
}

std::unique_ptr<DataSource> DataSourceBuilder::assemble() {
  // TODO (Michael Lettrich): Think of a nicer way to do assembly
  std::unique_ptr<SampleProvider> sampleProvider;
  std::unique_ptr<FileSampleProvider> fileSampleProvider;
  switch (fileType) {
    case ARFF:
      fileSampleProvider = std::make_unique<ArffFileSampleProvider>();
      break;
    case NONE:
      break;
    default:
      base::data_exception("Unknown file type");
      break;
  }

  if (isCompressed) {
    fileSampleProvider = std::make_unique<GzipFileSampleDecorator>(std::move(fileSampleProvider));
  }

  auto state = std::make_shared<DataSourceState>(filePath, batchSize, numBatches);
  sampleProvider = std::move(fileSampleProvider);
  return std::make_unique<DataSource>(state, std::move(sampleProvider));
}

} /* namespace datadriven */
} /* namespace sgpp */
