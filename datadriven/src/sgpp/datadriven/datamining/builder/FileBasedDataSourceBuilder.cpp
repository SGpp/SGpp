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

#include "FileBasedDataSourceBuilder.hpp"

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceState.hpp>
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

FileBasedDataSourceBuilder::FileBasedDataSourceBuilder()
    : DataSourceBuilder(), fileType(NONE), filePath(""), isCompressed(false) {}

FileBasedDataSourceBuilder::~FileBasedDataSourceBuilder() {}

DataSourceBuilder& FileBasedDataSourceBuilder::withFileType(std::string fileType) {
  if (arff.compare(fileType) == 0) {
    this->fileType = ARFF;
  } else {
    base::data_exception("Unknown file type");
  }
  return *this;
}

DataSourceBuilder& FileBasedDataSourceBuilder::withCompression(bool isCompressed) {
  this->isCompressed = isCompressed;
  return *this;
}

DataSourceBuilder& FileBasedDataSourceBuilder::withPath(std::string filePath) {
  this->filePath = filePath;
  grabTypeInfoFromFilePath();
  return *this;
}

std::unique_ptr<DataSource> FileBasedDataSourceBuilder::assemble() {
  std::unique_ptr<FileSampleProvider> fileSampleProvider;
  switch (fileType) {
    case ARFF:
      fileSampleProvider = std::make_unique<ArffFileSampleProvider>();
      break;
    case NONE:
    default:
      base::data_exception("Unknown file type");
      break;
  }

  if (isCompressed) {
    fileSampleProvider = std::make_unique<GzipFileSampleDecorator>(std::move(fileSampleProvider));
  }

  auto state = std::make_shared<DataSourceState>(filePath, numBatches, batchSize);
  std::unique_ptr<SampleProvider> sampleProvider = std::move(fileSampleProvider);
  return std::make_unique<DataSource>(state, std::move(sampleProvider));
}

void FileBasedDataSourceBuilder::grabTypeInfoFromFilePath() {
  // tokenize string
  std::vector<std::string> tokens;

  // split the string
  tokenize(filePath, ".", tokens);
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

// TODO(Michael Lettrich): move to own class
void FileBasedDataSourceBuilder::tokenize(const std::string& s, const char* delim,
                                          std::vector<std::string>& v) {
  // to avoid modifying original string
  // first duplicate the original string and return a char pointer then free the memory
  char* dup = strdup(s.c_str());
  char* token = strtok(dup, delim);
  while (token != NULL) {
    v.push_back(std::string(token));
    // the call is treated as a subsequent calls to strtok:
    // the function continues from where it left in previous invocation
    token = strtok(NULL, delim);
  }
  free(dup);
}

} /* namespace datadriven */
} /* namespace sgpp */
