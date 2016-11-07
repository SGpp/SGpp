// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NOTAKNOT_BSPLINE_MODIFIED_BASE_HPP
#define NOTAKNOT_BSPLINE_MODIFIED_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/NotAKnotBsplineBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * Not-a-knot B-spline basis.
 */
template <class LT, class IT>
class NotAKnotBsplineModifiedBasis: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NotAKnotBsplineModifiedBasis() : notAKnotBsplineBasis(NotAKnotBsplineBasis<LT, IT>()) {
  }

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NotAKnotBsplineModifiedBasis(size_t degree) :
        notAKnotBsplineBasis(NotAKnotBsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 3) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NotAKnotBsplineModifiedBasis() override {
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    double t = x * static_cast<double>(hInv) - static_cast<double>(i);
    
    switch (getDegree()) {
      case 1:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return std::max(1.0 - std::abs(t), 0.0);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          // l >= 3, i = 1
          return std::max(1.0 - t, 0.0);
        }

      case 3:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return notAKnotBsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l == 2, i == 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 1.0/40.0;
              result *= t;
              result = -3.0/5.0 + result * t;
              result = 6.0/5.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0/40.0;
              result = 3.0/20.0 + result * t;
              result = -3.0/10.0 + result * t;
              result = 1.0/5.0 + result * t;
              return result;
            }
          } else {
            // l >= 3, i == 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 1.0/24.0;
              result *= t;
              result = -3.0/4.0 + result * t;
              result = 5.0/4.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0/12.0;
              result = 1.0/4.0 + result * t;
              result = -1.0/4.0 + result * t;
              result = 1.0/12.0 + result * t;
              return result;
            }
          }
        }

      default:
        return 0.0;
    }
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const {
    return notAKnotBsplineBasis.getDegree();
  }

 protected:
  /// B-spline basis for B-spline evaluation
  NotAKnotBsplineBasis<LT, IT> notAKnotBsplineBasis;
};

// default type-def (unsigned int for level and index)
typedef NotAKnotBsplineModifiedBasis<unsigned int, unsigned int> SNotAKnotBsplineModifiedBase;

}  // namespace base
}  // namespace sgpp

#endif /* NOTAKNOT_BSPLINE_MODIFIED_BASE_HPP */
