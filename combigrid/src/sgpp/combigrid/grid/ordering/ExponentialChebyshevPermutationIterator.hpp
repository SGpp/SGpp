// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/common/AbstractPermutationIterator.hpp>

namespace sgpp {
namespace combigrid {

class ExponentialChebyshevPermutationIterator : public AbstractPermutationIterator {
  size_t currentIndex;
  size_t level;
  size_t numPoints;

 public:
  ExponentialChebyshevPermutationIterator(size_t level, size_t numPoints, size_t currentIndex = 0);

  virtual ~ExponentialChebyshevPermutationIterator();

  /**
   * Sets the iterator back to the start
   */
  virtual void reset();

  virtual size_t value();

  virtual void moveToNext();

  virtual std::shared_ptr<AbstractPermutationIterator> clone();
};

} /* namespace combigrid */
} /* namespace sgpp*/
