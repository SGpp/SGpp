// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/configuration/CrossvalidationConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctor.hpp>

#include <random>

namespace sgpp {
namespace datadriven {

/**
 * A pseudo shuffling functor, that maps indices 0...foldSize-1 to the current fold while mapping
 * the other indices sequentially to the rest of the dataset. Another shuffling functor is chained
 * to the results of this functor in order to retrieve a real shuffling.
 */
class DataShufflingFunctorCrossValidation : public DataShufflingFunctor {
 public:
  /**
   * Constructor.
   * @param crossValidationConfig configuration for the cross validation
   * @param shuffling the shuffling functor this functor is chained after
   */
  DataShufflingFunctorCrossValidation(const CrossvalidationConfiguration crossValidationConfig,
      DataShufflingFunctor* shuffling);

  /**
   * Clone pattern.
   * @return identical copy of this instance
   */
  DataShufflingFunctor* clone() const override;

  /**
   * Set the index of the current fold
   * @param fold the index of the current fold
   */
  void setFold(size_t fold);

  /**
   * Returns the size of the fold currently used for validation
   * @param numSamples the number of samples in total
   * @return size of current fold
   */
  size_t getCurrentFoldSize(size_t numSamples);

  /**
   * Overload the function-call operator that maps indexes to indexes via a permutation
   * of the entire index set. The permutation used is the identity.
   * @param idx the original index
   * @param numSamples the total number of indexes to permute
   * @return idx the index after the permutation (simply the input)
   */
  size_t operator()(size_t idx, size_t numSamples) override;

 private:
  /**
   * The shuffling functor this one is chained after
   */
  DataShufflingFunctor* shuffling;

  /**
   * Configuration for the cross validation
   */
  const CrossvalidationConfiguration crossValidationConfig;

  /**
   * Current fold
   */
  size_t currentFold;
};

} /* namespace datadriven */
} /* namespace sgpp */


