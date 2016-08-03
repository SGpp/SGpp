// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_ABSTRACTPERMUTATIONITERATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_ABSTRACTPERMUTATIONITERATOR_HPP_

#include <sgpp/globaldef.hpp>
#include <cstddef>
#include <memory>

namespace sgpp {
namespace combigrid {

class AbstractPermutationIterator {
 public:
  virtual ~AbstractPermutationIterator();

  /**
   * Sets the iterator back to the start
   */
  virtual void reset() = 0;

  virtual size_t value() = 0;

  virtual void moveToNext() = 0;

  virtual std::shared_ptr<AbstractPermutationIterator> clone() = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_ABSTRACTPERMUTATIONITERATOR_HPP_ */
