/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SampleProviderModule.hpp
 *
 *  Created on: 17.05.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceIterator.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/SampleProvider.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/globaldef.hpp>

#include <memory>
#include <string>

namespace sgpp {
namespace datadriven {

class DataSourceIterator;

class DataSource {
 public:
  DataSource(DataSourceConfig config, SampleProvider* sampleProvider);
  virtual ~DataSource();

  DataSourceConfig& getConfig();
  Dataset* getNextSamples();
  DataSourceIterator begin();
  DataSourceIterator end();
  size_t getCurrentIteration();

 protected:
  DataSourceConfig config;
  size_t currentIteration;
  std::shared_ptr<SampleProvider> sampleProvider;
};

} /* namespace datadriven */
} /* namespace sgpp */
