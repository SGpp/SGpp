/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SampleProvider.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun, Michael Lettrich
 */

#pragma once

#include <memory>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

class SampleProvider {
 public:
  SampleProvider(){};
  virtual ~SampleProvider(){};

  /**
   * Selects a certain number of samples
   * @param howMany number of samples to return
   * @return Dataset* pointer to a Dataset containing all samples.
   */
  virtual Dataset* getNextSamples(size_t howMany) = 0;

  /**
   * Returns all samples
   * @return std::unique_ptr<Dataset> A smart pointer to a Dataset.
   */
  virtual Dataset* getAllSamples() = 0;

  /**
   * Returns the dimensionality of the data source
   * @return dimensionality of the dataset.
   */
  virtual size_t getDim() = 0;
};
}
}
