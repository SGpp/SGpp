// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctor.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Concrete factory to build an instance of #sgpp::datadriven::DataShufflingFunctor
 */
class DataShufflingFunctorFactory {
 public:
  /**
   * Default constructor
   */
  DataShufflingFunctorFactory() = default;

  /**
   * Create an instance of a #sgpp::datadriven::DataShufflingFunctor object based on the
   * configuration
   * @param config configuration for the data source
   * @return Fully configured instance of a  #sgpp::datadriven::DataShufflingFunctor object.
   */
  DataShufflingFunctor* buildDataShufflingFunctor(const DataSourceConfig& config) const;
};

} /* namespace datadriven */
} /* namespace sgpp */
