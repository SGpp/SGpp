/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FileSampleProviderBuilder.hpp
 *
 *  Created on: 16.05.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>

#include <memory>
#include <string>

enum FileType { NONE, ARFF };

namespace sgpp {
namespace datadriven {
// TODO(Michael Lettrich): parse filetype?
class DataSourceBuilder {
 public:
  DataSourceBuilder();
  virtual ~DataSourceBuilder();
  DataSourceBuilder& inBatches(size_t howMany);
  DataSourceBuilder& withBatchSize(size_t batchSize);
  virtual std::unique_ptr<DataSource> assemble() = 0;

 protected:
  size_t batchSize;
  size_t numBatches;
};

} /* namespace datadriven */
} /* namespace sgpp */
