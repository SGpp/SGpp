// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEAR_BASE_HPP
#define LINEAR_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <algorithm>

namespace sgpp {
namespace base {

/**
 * Linear basis on Noboundary grids.
 */
template <class LT, class IT>
class LinearBasis : public Basis<LT, IT> {
 public:
  /**
   * Destructor.
   */
  ~LinearBasis() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of linear basis function
   */
  inline double eval(LT l, IT i, double x) override {
    return std::max(
        1.0 - std::abs(static_cast<double>(static_cast<IT>(1) << l) * x - static_cast<double>(i)),
        0.0);
  }

  size_t getDegree() const override { return 1; }
};

// default type-def (unsigned int for level and index)
typedef LinearBasis<unsigned int, unsigned int> SLinearBase;

}  // namespace base
}  // namespace sgpp

#endif /* LINEAR_BASE_HPP */
