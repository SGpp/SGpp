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

#include "DataSourceBuilder.hpp"

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>

namespace sgpp {
namespace datadriven {

class FileBasedDataSourceBuilder : public DataSourceBuilder {
 public:
  FileBasedDataSourceBuilder();
  virtual ~FileBasedDataSourceBuilder();
  DataSourceBuilder& withFileType(std::string fileType);
  DataSourceBuilder& withCompression(bool isCompressed);
  DataSourceBuilder& withPath(std::string filePath);
  virtual std::unique_ptr<DataSource> assemble();

 protected:
  FileType fileType;
  std::string filePath;
  bool isCompressed;

 private:
  void grabTypeInfoFromFilePath();
  void tokenize(const std::string& s, const char* delim, std::vector<std::string>& v);
};

} /* namespace datadriven */
} /* namespace sgpp */
