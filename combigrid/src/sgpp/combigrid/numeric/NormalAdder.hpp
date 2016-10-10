// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NORMALADDER_HPP_
#define NORMALADDER_HPP_
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace combigrid {

/**
 * Standard adder with the same interface as the Kahan adder.
 */
class NormalAdder {
  double sum = 0.0;

 public:
  void add(double x) { sum += x; }

  double value() const { return sum; }
};

}  // namespace combigrid
} /* namespace sgpp*/

#endif /* NORMALADDER_HPP_ */
