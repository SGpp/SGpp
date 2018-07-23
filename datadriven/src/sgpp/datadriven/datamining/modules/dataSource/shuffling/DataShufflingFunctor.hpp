/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataShufflingFunctor.hpp
 *
 *  Created on: Jul 20, 2018
 *      Author: dominik
 */

#pragma once

#include <sgpp/globaldef.hpp>

#include <random>

namespace sgpp {
namespace datadriven {

/**
 * A class to provide functionality to shuffle (reorder) the data samples before the sample
 * provider accesses it. This is neccessary for cross validation and general shuffling of data.
 */
class DataShufflingFunctor {
 public:
  DataShufflingFunctor();
  virtual ~DataShufflingFunctor() = default;
  /**
   * Overload the function-call operator that maps indexes to indexes via a permutation
   * of the entire index set.
   * @param idx the original index
   * @return idx the index after the permutation
   */
  virtual size_t operator()(size_t idx);
};

} /* namespace datadriven */
} /* namespace sgpp */

