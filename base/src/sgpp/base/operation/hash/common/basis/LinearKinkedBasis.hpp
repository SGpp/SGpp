// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>

namespace sgpp {
namespace base {

/**
 * Kinked linear basis on Noboundary grids.
 */
template <class LT, class IT>
class LinearKinkedBasis : public Basis<LT, IT> {
 public:
  /**
   * Destructor.
   */
  ~LinearKinkedBasis() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of kinked linear basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);

    if (l == 1) {
      // first level
      return 1.0;
    } else if (i == 1) {
      // left kinked basis function
      return ((x <= 2.0 / hInvDbl) ? std::min((2.0 - hInvDbl * x), 1.0) : 0.0);
    } else if (i == hInv - 1) {
      // right kinked basis function
      return ((x >= 1.0 - 2.0 / hInvDbl) ?
      std::min((hInvDbl * x - static_cast<double>(i) + 1.0), 1.0) : 0.0);
    } else {
      // interior basis function
      return std::max(1.0 - std::abs(hInvDbl * x - static_cast<double>(i)), 0.0);
    }
  }

  inline size_t getDegree() const override { return 1; }

  double getIntegral(LT level, IT index) override {
    const IT hInv = static_cast<IT>(1) << level;

    if (level == 1) {
      // first level
      return 1.0;
    } else if ((index == 1) || (index == hInv - 1)) {
      // left and right kinked basis functions
      return 1.5 / hInv;
    } else {
      // interior basis function
      return 1. / hInv;
    }
  }
};

// default type-def (unsigned int for level and index)
typedef LinearKinkedBasis<unsigned int, unsigned int> SLinearKinkedBase;

}  // namespace base
}  // namespace sgpp
