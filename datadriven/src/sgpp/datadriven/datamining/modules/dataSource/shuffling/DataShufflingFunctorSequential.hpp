// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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
class DataShufflingFunctorSequential : public DataShufflingFunctor {
 public:
  /**
   * Clone pattern.
   * @return identical copy of this instance
   */
  DataShufflingFunctor* clone() const override;

  /**
   * Overload the function-call operator that maps indexes to indexes via a permutation
   * of the entire index set. The permutation used is the identity.
   * @param idx the original index
   * @param numSamples the total number of indexes to permute
   * @return idx the index after the permutation (simply the input)
   */
  size_t operator()(size_t idx, size_t numSamples) override;
};

} /* namespace datadriven */
} /* namespace sgpp */

