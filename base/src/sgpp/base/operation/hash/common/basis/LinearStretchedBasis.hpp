// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARSTRETCHED_BASE_HPP
#define LINEARSTRETCHED_BASE_HPP

#include <sgpp/base/grid/common/Stretching.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <algorithm>


namespace sgpp {
namespace base {

/**
 * linearstretched base functions.
 * And here we have another implicit dependence on tensor products
 */
template<class LT, class IT>
class LinearStretchedBasis: public LinearBasis<LT, IT> {
 public:
  /**
   * Destructor.
   */
  ~LinearStretchedBasis() override {
  }

  /**
   * Evaluate a basis function.
   * Has a dependence on the absolute position of grid point and support.

  double eval(LT level, IT index, double p)
  {
    return 1.0 - fabs((1<<level) * p - index);
  }

   *
   * Evaluate a basis function.
   * Has a dependence on the absolute position of grid point and support.
   *
   * This version catches errors, that occur if a basis function
   * is evaluated outside its domain

  double evalSave(LT level, IT index, double p)
  {
    return std::max(1.0 - fabs((1<<level) * p - index), 0.0);
  }*/

  /*
   * evaluate a basis function
   * Has a dependence on the position of two grid points with values 1 and 0 and the
   * support position
   */
  double stretchedEval(double p, double pos0, double pos1) {
    return (p - pos0) / (pos1 - pos0);
  }
};

// default type-def (unsigned int for level and index)
typedef LinearStretchedBasis<unsigned int, unsigned int> SLinearStretchedBase;

}  // namespace base
}  // namespace sgpp

#endif /* LINEARSTRETCHED_BASE_HPP */
