// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARPERIODICBASE_HPP
#define LINEARPERIODICBASE_HPP

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>

namespace sgpp {
namespace base {

/**
 * linear basis functions with boundaries
 * And here we have another implicit dependence on tensor products
 *
 */
template <class LT, class IT>
class LinearPeriodicBasis : public Basis<LT, IT> {
 public:
  /**
   * Destructor.
   */
  ~LinearPeriodicBasis() override {}

  /**
   * Evaluate a basis function.
   * Has a dependence on the absolute position of grid point and support.
   *
   * @param l the level of the current basis function
   * @param i the index of the current basis function
   * @param x the absolute position of the evaluation point
   */
  inline double eval(LT l, IT i, double x) override {
    if (l == 0) {
      return fabs(2 * x - 1);
    } else {
      return std::max(1.0 - fabs(static_cast<double>(1 << l) * x - static_cast<double>(i)), 0.0);
    }
  }

  /**
   * Evaluate a basis function with an offset and scaling factor
   * Has a dependence on the absolute position of grid point and support.
   *
   * @param level the level of the current basis function
   * @param index the index of the current basis function
   * @param p the absolute position of the evaluation point
   * @param q the scaling factor of the basis function
   * @param t the offset of the basis function
   */
  double eval(LT level, IT index, double p, double q, double t) {
    if (level == 0) {
      return ((1.0 / q) * (fabs((2 * (p - t)) - (q))));
    } else {
      return std::max(1.0 - ((1.0 / q) * (fabs((static_cast<double>(1 << level) * (p - t)) -
                                               (q * static_cast<double>(index))))),
                      0.0);
    }
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "LinearPeriodicBasis: evalDx not implemented" << std::endl;
    return -1;
  }

  inline double getIntegral(LT level, IT index) override { return -1.0; }

  inline size_t getDegree() const override { return 1; }
};

// default type-def (unsigned int for level and index)
typedef LinearPeriodicBasis<unsigned int, unsigned int> SLinearPeriodicBasis;

}  // namespace base
}  // namespace sgpp

#endif /* LINEARPERIODICBASE_HPP */
