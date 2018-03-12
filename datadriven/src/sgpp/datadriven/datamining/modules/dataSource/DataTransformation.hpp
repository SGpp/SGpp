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

#ifndef DATATRANSFORMATION_HPP
#define DATATRANSFORMATION_HPP

#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>

namespace sgpp {
namespace datadriven {

/**
 * DataTransformation is an abstraction for an object that provides different transformations on
 * datasets, for example Rosenblatt-transformation to get a uniform distribution over the unit cube.
 */

class DataTransformation {
 public:
  /**
   * Default constructor
   */
  DataTransformation() = default;

  /**
     * Virtual destructor
     */
  virtual ~DataTransformation() = default;

  /**
   * Performs a data transformation on a given dataset for a data transformationn built with
   * DataTransformationBuilder
   *
   * @param dataset pointer to the dataset to be transformed
   * @return pointer to the transformed dataset
   */
  virtual Dataset* doTransformation(Dataset* dataset) = 0;

  /**
     * Performs the backwards transformation on a given dataset for a data transformationn
     * built with DataTransformationBuilder
     *
     * @param dataset pointer to the dataset to be transformed backwards
     * @return pointer to the backwards transformed dataset
     */
  virtual Dataset* doInverseTransformation(Dataset* dataset) = 0;
};
} /* namespace datadriven */
} /* namespace sgpp */
#endif /* DATATRANSFORMATION_HPP */
