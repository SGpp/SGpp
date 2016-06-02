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

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceIterator.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceState.hpp>
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
  DataSource(std::shared_ptr<DataSourceState> state,
             std::unique_ptr<SampleProvider> sampleProvider);
  DataSource(DataSource&& ds);
  virtual ~DataSource();

  std::unique_ptr<Dataset> getNextSamples();
  DataSourceIterator begin();
  DataSourceIterator end();

 protected:
  std::shared_ptr<DataSourceState> state;
  std::unique_ptr<SampleProvider> sampleProvider;
};

} /* namespace datadriven */
} /* namespace sgpp */
