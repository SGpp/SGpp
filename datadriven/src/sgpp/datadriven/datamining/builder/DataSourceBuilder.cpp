/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FileBasedDataSourceBuilder.cpp
 *
 *  Created on: 01.06.2016
 *      Author: Michael Lettrich
 */

#include "DataSourceBuilder.hpp"

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/base/StringTokenizer.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/GzipFileSampleDecorator.hpp>

#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

const std::string arff = "arff";
const std::string gz = "gz";

namespace sgpp {
namespace datadriven {

DataSourceBuilder::DataSourceBuilder() : config(), fileType(NONE), isCompressed(false) {
  config.numBatches = 1;
}

DataSourceBuilder::~DataSourceBuilder() {}

DataSourceBuilder& DataSourceBuilder::withFileType(std::string fileType) {
  if (arff.compare(fileType) == 0) {
    this->fileType = ARFF;
  } else {
    base::data_exception("Unknown file type");
  }
  return *this;
}

DataSourceBuilder& DataSourceBuilder::inBatches(size_t howMany) {
  config.numBatches = howMany;
  return *this;
}

DataSourceBuilder& DataSourceBuilder::withBatchSize(size_t batchSize) {
  config.batchSize = batchSize;
  return *this;
}

DataSourceBuilder& DataSourceBuilder::withCompression(bool isCompressed) {
  this->isCompressed = isCompressed;
  return *this;
}

DataSourceBuilder& DataSourceBuilder::withPath(std::string filePath) {
  config.filePath = filePath;
  grabTypeInfoFromFilePath();
  return *this;
}

DataSource* DataSourceBuilder::assemble() {
  SampleProvider* sampleProvider;

  switch (fileType) {
    case ARFF:
      sampleProvider = new ArffFileSampleProvider;
      break;
    case NONE:
    default:
      base::data_exception("Unknown file type");
      break;
  }

  if (isCompressed) {
    sampleProvider = new GzipFileSampleDecorator(static_cast<FileSampleProvider*>(sampleProvider));
  }

  return new DataSource(config, sampleProvider);
}

void DataSourceBuilder::grabTypeInfoFromFilePath() {
  // tokenize string
  std::vector<std::string> tokens;

  // split the string
  StringTokenizer::tokenize(config.filePath, ".", tokens);
  // convert to lower case
  for (auto t : tokens) {
    // TODO(Michael Lettrich): test if this works with umlauts
    std::transform(t.begin(), t.end(), t.begin(), ::tolower);
  }
  // check if there is gz compression
  if (tokens.back().compare(gz) == 0) {
    isCompressed = true;
  }

  // check if we can find ARFF
  if (std::find(tokens.begin(), tokens.end(), arff) != tokens.end()) {
    fileType = ARFF;
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
