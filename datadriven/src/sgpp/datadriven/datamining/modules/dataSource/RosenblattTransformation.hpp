/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * RosenblattTransformation.hpp
 *
 *  Created on: 22.01.2018
 *      Author: Lars Wolfsteller
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformation.hpp>
#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

namespace sgpp {

using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::DataVector;
using sgpp::base::DataMatrix;
using sgpp::optimization::RandomNumberGenerator;

namespace datadriven {

class RosenblattTransformation : public DataTransformation {
 public:
  /**
   * Default constructor
   */
  RosenblattTransformation(Dataset* dataset, size_t numSamples);

  Dataset* doTransformation(Dataset* dataset);
  Dataset* doInverseTransformation(Dataset* dataset);

 private:
  /**
   * the sparse grid that approximates the data.
   */
  sgpp::base::Grid* grid;

  /**
   * hierarchical surpluses of the #grid.
   */
  sgpp::base::DataVector* alpha;

  /**
   * Pointer to #sgpp::datadriven::Dataset
   */
  Dataset* datasetTransformed;

  /**
   * Pointer to #sgpp::datadriven::Dataset
   */
  Dataset* datasetInvTransformed;

  /**
   * Number of samples for calculation of pdf / alpha, default: 1000
   */
  size_t numSamples = 1000;
};
} /* namespace datadriven */
} /* namespace sgpp */
