// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#ifndef DATATRANSFORMATIONBUILDER_HPP
#define DATATRANSFORMATIONBUILDER_HPP

#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformation.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Helper class to build all kinds of transformation based on given configuration
 */

class DataTransformationBuilder {
 public:
  DataTransformationBuilder() = default;
  virtual ~DataTransformationBuilder() = default;

  /**
   * Initializes a data transformation for a given dataset according to the
   * specified configuration
   *
   * @param config DataSourceConfig containing DataTransformationType and
   *  numSamplesForTranformation to calculate probability density function pdf / alpha
   * @return pointer to the initialized dataTransformation
   */
  DataTransformation *buildTransformation(DataTransformationConfig config);
};
} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATATRANSFORMATIONBUILDER_HPP */
