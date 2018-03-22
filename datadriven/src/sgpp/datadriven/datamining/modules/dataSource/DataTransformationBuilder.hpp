/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataTransformation.hpp
 *
 *  Created on: 22.01.2018
 *      Author: Lars Wolfsteller
 */

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
   * @param dataset pointer to the dataset on which transformation will be performed
   * @return pointer to the initialized dataTransformation
   */
  DataTransformation* buildTransformation(DataTransformationConfig config, Dataset* dataset);
};
} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATATRANSFORMATIONBUILDER_HPP */
