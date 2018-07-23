/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataShufflingFunctorSequential.hpp
 *
 *  Created on: Jul 20, 2018
 *      Author: dominik
 */

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctor.hpp>

#include <random>

namespace sgpp {
namespace datadriven {

/**
 * A simple shuffling functor for data samples that performs no shuffling at all, i.e. the
 * permutation used is the identity itself.
 */
class DataShufflingFunctorSequential : DataShufflingFunctor {
 public:
  /**
   * Overload the function-call operator that maps indexes to indexes via a permutation
   * of the entire index set. The permutation used is the identity.
   * @param idx the original index
   * @return idx the index after the permutation (simply the input)
   */
  size_t operator()(size_t idx) override;
};

} /* namespace datadriven */
} /* namespace sgpp */

