// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NOTAKNOT_BSPLINE_BASE_HPP
#define NOTAKNOT_BSPLINE_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * Not-a-knot B-spline basis.
 */
template <class LT, class IT>
class NotAKnotBsplineBasis: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NotAKnotBsplineBasis() : bsplineBasis(BsplineBasis<LT, IT>()) {
  }

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NotAKnotBsplineBasis(size_t degree) :
        bsplineBasis(BsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 3) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NotAKnotBsplineBasis() override {
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
        return std::max(1.0 - std::abs(t), 0.0);

      case 3:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return 1.0 - x;
          } else {
            // l = 0, i = 1
            return x;
          }
        } else if (l == 1) {
          // l = 1, i = 1
          return 1.0 - t * t;
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 4, 3 < i < 2^l - 3
          return bsplineBasis.eval(l, i, x);
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
              double result = 1.0/10.0;
              result = -9.0/20.0 + result * t;
              result = 3.0/10.0 + result * t;
              result = 3.0/5.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0/40.0;
              result = 3.0/20.0 + result * t;
              result = -3.0/10.0 + result * t;
              result = 1.0/5.0 + result * t;
              return result;
            }
          } else if (i == 1) {
            // l >= 3, i == 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 1.0/8.0;
              result = -1.0/2.0 + result * t;
              result = 1.0/4.0 + result * t;
              result = 7.0/12.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0/12.0;
              result = 1.0/4.0 + result * t;
              result = -1.0/4.0 + result * t;
              result = 1.0/12.0 + result * t;
              return result;
            }
          } else {
            // l >= 3, i == 3
            if ((t < -3.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 3.0;
              double result = 1.0/24.0;
              result *= t;
              result *= t;
              result *= t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -3.0/8.0;
              result = 1.0/4.0 + result * t;
              result = 1.0/2.0 + result * t;
              result = 1.0/3.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 11.0/24.0;
              result = -7.0/8.0 + result * t;
              result = -1.0/8.0 + result * t;
              result = 17.0/24.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0/6.0;
              result = 1.0/2.0 + result * t;
              result = -1.0/2.0 + result * t;
              result = 1.0/6.0 + result * t;
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
    return bsplineBasis.getDegree();
  }

 protected:
  /// B-spline basis for B-spline evaluation
  BsplineBasis<LT, IT> bsplineBasis;
};

// default type-def (unsigned int for level and index)
typedef NotAKnotBsplineBasis<unsigned int, unsigned int> SNotAKnotBsplineBase;

}  // namespace base
}  // namespace sgpp

#endif /* NOTAKNOT_BSPLINE_BASE_HPP */
