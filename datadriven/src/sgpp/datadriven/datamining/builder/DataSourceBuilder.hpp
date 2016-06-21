/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FileBasedDataSourceBuilder.hpp
 *
 *  Created on: 01.06.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>

#include <memory.h>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

enum FileType { NONE, ARFF };

class DataSourceBuilder {
 public:
  DataSourceBuilder();
  virtual ~DataSourceBuilder();
  DataSourceBuilder& withPath(std::string filePath);
  DataSourceBuilder& withCompression(bool isCompressed);
  DataSourceBuilder& withFileType(std::string fileType);
  DataSourceBuilder& inBatches(size_t howMany);
  DataSourceBuilder& withBatchSize(size_t batchSize);
  virtual std::unique_ptr<DataSource> assemble();

 private:
  void grabTypeInfoFromFilePath();
  void tokenize(const std::string& s, const char* delim, std::vector<std::string>& v);

  FileType fileType;
  std::string filePath;
  bool isCompressed;
  size_t batchSize;
  size_t numBatches;
};

} /* namespace datadriven */
} /* namespace sgpp */
