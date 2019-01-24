// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEAR_MODIFIED_BASE_HPP
#define LINEAR_MODIFIED_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>

namespace sgpp {
namespace base {

/**
 * Modified linear basis on Noboundary grids.
 */
template <class LT, class IT>
class LinearModifiedBasis : public Basis<LT, IT> {
 public:
  /**
   * Destructor.
   */
  ~LinearModifiedBasis() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of modified linear basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);

    if (l == 1) {
      // first level
      return 1.0;
    } else if (i == 1) {
      // left modified basis function
      return ((x <= 2.0 / hInvDbl) ? (2.0 - hInvDbl * x) : 0.0);
    } else if (i == hInv - 1) {
      // right modified basis function
      return ((x >= 1.0 - 2.0 / hInvDbl) ? (hInvDbl * x - static_cast<double>(i) + 1.0) : 0.0);
    } else {
      // interior basis function
      return std::max(1.0 - std::abs(hInvDbl * x - static_cast<double>(i)), 0.0);
    }
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "LinearModifiedBasis: evalDx not implemented" << std::endl;
    return -1;
  }

  inline size_t getDegree() const override { return 1; }

  double getIntegral(LT level, IT index) override {
    const IT hInv = static_cast<IT>(1) << level;

    if (level == 1) {
      // first level
      return 1.0;
    } else if ((index == 1) || (index == hInv - 1)) {
      // left and right modified basis functions
      return 2. / hInv;
    } else {
      // interior basis function
      return 1. / hInv;
    }
  }
};

// default type-def (unsigned int for level and index)
typedef LinearModifiedBasis<unsigned int, unsigned int> SLinearModifiedBase;

}  // namespace base
}  // namespace sgpp

#endif /* LINEAR_MODIFIED_BASE_HPP */
