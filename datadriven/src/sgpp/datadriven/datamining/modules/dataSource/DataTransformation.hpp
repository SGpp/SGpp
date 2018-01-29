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

#include <sgpp/datadriven/tools/Dataset.hpp>

namespace sgpp {
namespace datadriven {

/**
 * DataTransformation is an abstraction for an object that provides different transformations on
 * datasets, for example Rosenblatt-transformation to get a uniform distribution over the unit cube.
 */

class DataTransformation {
 public:
  DataTransformation();
  virtual Dataset* doTransformation(Dataset* dataset);
  virtual Dataset* doInverseTransformation(Dataset* dataset);
  virtual ~DataTransformation();
};
} /* namespace datadriven */
} /* namespace sgpp */
