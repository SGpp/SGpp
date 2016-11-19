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

#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/globaldef.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

/**
 * SampleProvider provides a unified interface for all objects that can supply #Dataset objects.
 * Datasets can be generated out of files (e.g. by reading ARFF or CSV files), functions
 * generating artificial datasets (e.g. Friedmann functions) or streams (e.g. from an open TCP/IP
 * connection) to name just a few.
 *
 * A Sample Provider is only responsible for getting data.
 */
class SampleProvider {
 public:
  SampleProvider(){};
  virtual ~SampleProvider(){};

  /**
   * Request a #Dataset containing the desired amount of samples. This interface is of primary
   * interest for streaming algorithms and batch learning, where datasets do not fit into memory and
   * thus have to be split in parts or functions generating artificial #Datasets.
   * @param howMany amount of samples requested. If callee can not fulfill the request, less points
   * are returned. It is the callers responsibility to check the size.
   * @return std::unique_ptr<Dataset> unique pointer to a #Dataset containing the samples, owned by
   * caller.
   */
  virtual auto getNextSamples(size_t howMany) -> std::unique_ptr<Dataset> = 0;

  /**
   * Return  a dataset containing all samples the #SampleProvider has to offer. This interface is
   * meant for reading entire datasets out of files.
   * @return std::unique_ptr<Dataset> unique pointer to a #Dataset containing the samples, owned by
   * caller.
   */
  virtual auto getAllSamples() -> std::unique_ptr<Dataset> = 0;
};
}
}
