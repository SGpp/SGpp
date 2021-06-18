// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEAR_BOUNDARY_BASE_HPP
#define LINEAR_BOUNDARY_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>

namespace sgpp {
namespace base {

/**
 * Linear basis on Boundary grids.
 */
template <class LT, class IT>
class LinearBoundaryBasis : public Basis<LT, IT> {
 public:
  /**
   * Destructor.
   */
  ~LinearBoundaryBasis() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of boundary linear basis function
   */
  inline double eval(LT l, IT i, double x) override {
    if (l == 0) {
      // first level
      if (i == 0) {
        return 1.0 - x;
      } else {
        return x;
      }
    } else {
      return std::max(
          1.0 - std::abs(static_cast<double>(static_cast<IT>(1) << l) * x - static_cast<double>(i)),
          0.0);
    }
  }

  /**
   * Evaluate basis function with offset and scaling factor.
   *
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @param q     scaling factor of basis function
   * @param t     offset of basis function
   */
  inline virtual double eval(LT l, IT i, double x, double q, double t) {
    return eval(l, i, (x - t) / q);
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "LinearBoundaryBasis: evalDx not implemented" << std::endl;
    return -1;
  }

  double getIntegral(LT level, IT index) override {
    if (level == 0) {
      return 0.5;
    } else {
      return 1. / static_cast<double>(static_cast<IT>(1) << level);
    }
  }

  inline size_t getDegree() const override { return 1; }
};

// default type-def (unsigned int for level and index)
typedef LinearBoundaryBasis<unsigned int, unsigned int> SLinearBoundaryBase;

}  // namespace base
}  // namespace sgpp

#endif /* LINEAR_BOUNDARY_BASE_HPP */
