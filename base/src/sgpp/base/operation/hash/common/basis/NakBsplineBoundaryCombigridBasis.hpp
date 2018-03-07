// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>

namespace sgpp {
namespace base {

/**
 * Hierarchical Not-a-knot B-spline basis. This basis is designed to represent Bspline boundary
 * interpolants from the combigrid module.
 * Therefore it has the following unusual choice of basis functions:
 * linear and quadratic terms on level 0
 * constant term on level 1
 *
 * This means that this basis of level 0 cannot represent constant functions! It is therefore not
 * suitabe for any application that cannot guarantee the existence of level 1 in every dimension
 *
 * A combigrid interpolant which shall be interpolated with this class needs level (1,...,1)
 * otherwise the boundary points cannot match
 */

template <class LT, class IT>
class NakBsplineBoundaryCombigridBasis : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakBsplineBoundaryCombigridBasis() : bsplineBasis(BsplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineBoundaryCombigridBasis(size_t degree)
      : bsplineBasis(BsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("NakBsplineBoundaryCombigridBasis: Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakBsplineBoundaryCombigridBasis() override {}

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

      // degree 3: global polynomials in x on Level 0 and 1, nak Bsplines from Level 3 on
      case 3:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return 1 - x * x;
          } else {
            // l = 0, i = 1
            return x;
          }
        } else if (l == 1) {
          if (i == 0) {
            // l = 1, i = 0
            return 0.5 * x * x - 1.5 * x + 1.0;
          } else if (i == 1) {
            // l = 1, i = 1
            return 1;
          } else {
            // l = 1, i = 2
            return 0.5 * x * x + 1.5 * x + 1.0;
          }
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 4, 3 < i < 2^l - 3
          return bsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (i == 0) {
            // l >= 2, i = 0
            if ((t < 0.0) || (t > 2.0)) {
              return 0.0;
            } else {
              double result = -4.1666666666666664e-02;
              result = 2.5000000000000000e-01 + result * t;
              result = -5.0000000000000000e-01 + result * t;
              result = 3.3333333333333331e-01 + result * t;
              return result;
            }
          } else if ((l == 2) && (i == 1)) {
            // l = 2, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 1.0 / 10.0;
              result = -9.0 / 20.0 + result * t;
              result = 3.0 / 10.0 + result * t;
              result = 3.0 / 5.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0 / 40.0;
              result = 3.0 / 20.0 + result * t;
              result = -3.0 / 10.0 + result * t;
              result = 1.0 / 5.0 + result * t;
              return result;
            }

          } else if (l == 2) {
            // l = 2, i = 2
            if ((t < -2.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 2.0;
              double result = -8.3333333333333329e-02;
              result = 2.0000000000000001e-01 + result * t;
              result = 2.0000000000000001e-01 + result * t;
              result = 6.6666666666666666e-02 + result * t;
              return result;
            } else {
              double result = 8.3333333333333329e-02;
              result = -2.9999999999999999e-01 + result * t;
              result *= t;
              result = 5.9999999999999998e-01 + result * t;
              return result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 1.0 / 8.0;
              result = -1.0 / 2.0 + result * t;
              result = 1.0 / 4.0 + result * t;
              result = 7.0 / 12.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0 / 12.0;
              result = 1.0 / 4.0 + result * t;
              result = -1.0 / 4.0 + result * t;
              result = 1.0 / 12.0 + result * t;
              return result;
            }
          } else if (i == 2) {
            // l >= 3, i = 2
            if ((t < -2.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 2.0;
              double result = -1.2500000000000000e-01;
              result = 2.5000000000000000e-01 + result * t;
              result = 2.5000000000000000e-01 + result * t;
              result = 8.3333333333333329e-02 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 2.9166666666666669e-01;
              result = -5.0000000000000000e-01 + result * t;
              result = -2.5000000000000000e-01 + result * t;
              result = 5.8333333333333337e-01 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.2500000000000000e-01;
              result = 3.7500000000000000e-01 + result * t;
              result = -3.7500000000000000e-01 + result * t;
              result = 1.2500000000000000e-01 + result * t;
              return result;
            }
          } else {
            // l >= 3, i = 3
            if ((t < -3.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 3.0;
              double result = 1.0 / 24.0;
              result *= t;
              result *= t;
              result *= t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -3.0 / 8.0;
              result = 1.0 / 4.0 + result * t;
              result = 1.0 / 2.0 + result * t;
              result = 1.0 / 3.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 11.0 / 24.0;
              result = -7.0 / 8.0 + result * t;
              result = -1.0 / 8.0 + result * t;
              result = 17.0 / 24.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0 / 6.0;
              result = 1.0 / 2.0 + result * t;
              result = -1.0 / 2.0 + result * t;
              result = 1.0 / 6.0 + result * t;
              return result;
            }
          }
        }

      // degree 5: Levels 0,1 and 2 polynomials, nak Bsplines from Level 3 on
      case 5:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return 1 - x * x;
          } else {
            // l = 0, i = 1
            return x;
          }
        } else if (l == 1) {
          if (i == 1) {
            // l = 1, i = 1
            return 1;
          }
        } else if (l == 2) {
          if (i == 1) {
            // l = 2, i = 1 : cubic polynomial, 0 in 0,0.5,0.75 and 1 in 0.25
            return 32 * x * (x - 0.5) * (x - 0.75);
          } else if (i == 3) {
            // l = 2, i = 3 : quartic polynomial, 0 in 0,0.25,0.5,1 and 1 in 0.75
            return x * x * x * x;  // x * (x - 0.25) * (x - 0.5) * (x - 1) * (-128.0 / 3.0);
          }
        } else if ((i > 5) && (i < hInv - 5)) {
          // l >= 4, 5 < i < 2^l - 5
          return bsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 107.0 / 30240.0;
              result = -17.0 / 756.0 + result * t;
              result = 1.0 / 378.0 + result * t;
              result = 37.0 / 378.0 + result * t;
              result = 109.0 / 756.0 + result * t;
              result = 253.0 / 3780.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -397.0 / 30240.0;
              result = 185.0 / 6048.0 + result * t;
              result = 155.0 / 3024.0 + result * t;
              result = -415.0 / 3024.0 + result * t;
              result = -1165.0 / 6048.0 + result * t;
              result = 2965.0 / 6048.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 233.0 / 30240.0;
              result = -53.0 / 1512.0 + result * t;
              result = 8.0 / 189.0 + result * t;
              result = 13.0 / 189.0 + result * t;
              result = -97.0 / 378.0 + result * t;
              result = 433.0 / 1890.0 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.0 / 4320.0;
              result = 1.0 / 288.0 + result * t;
              result = -1.0 / 48.0 + result * t;
              result = 1.0 / 16.0 + result * t;
              result = -3.0 / 32.0 + result * t;
              result = 9.0 / 160.0 + result * t;
              return result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 1.0 / 504.0;
              result = -1.0 / 42.0 + result * t;
              result = 2.0 / 21.0 + result * t;
              result = -2.0 / 21.0 + result * t;
              result = -5.0 / 21.0 + result * t;
              result = 47.0 / 105.0 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.0 / 840.0;
              result = 1.0 / 168.0 + result * t;
              result = -1.0 / 84.0 + result * t;
              result = 1.0 / 84.0 + result * t;
              result = -1.0 / 168.0 + result * t;
              result = 1.0 / 840.0 + result * t;
              return result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 1.0 / 252.0;
              result = -1.0 / 42.0 + result * t;
              result *= t;
              result = 2.0 / 21.0 + result * t;
              result = 1.0 / 7.0 + result * t;
              result = 1.0 / 15.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -23.0 / 1260.0;
              result = 1.0 / 28.0 + result * t;
              result = 1.0 / 14.0 + result * t;
              result = -5.0 / 42.0 + result * t;
              result = -1.0 / 4.0 + result * t;
              result = 163.0 / 420.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 19.0 / 1260.0;
              result = -1.0 / 18.0 + result * t;
              result = 2.0 / 63.0 + result * t;
              result = 8.0 / 63.0 + result * t;
              result = -2.0 / 9.0 + result * t;
              result = 34.0 / 315.0 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.0 / 252.0;
              result = 5.0 / 252.0 + result * t;
              result = -5.0 / 126.0 + result * t;
              result = 5.0 / 126.0 + result * t;
              result = -5.0 / 252.0 + result * t;
              result = 1.0 / 252.0 + result * t;
              return result;
            }
          } else {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < -2.0) {
              t += 5.0;
              double result = 1.0 / 2520.0;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = -11.0 / 504.0;
              result = 1.0 / 168.0 + result * t;
              result = 1.0 / 28.0 + result * t;
              result = 3.0 / 28.0 + result * t;
              result = 9.0 / 56.0 + result * t;
              result = 27.0 / 280.0 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 31.0 / 504.0;
              result = -13.0 / 126.0 + result * t;
              result = -10.0 / 63.0 + result * t;
              result = 2.0 / 63.0 + result * t;
              result = 25.0 / 63.0 + result * t;
              result = 121.0 / 315.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -181.0 / 2520.0;
              result = 103.0 / 504.0 + result * t;
              result = 11.0 / 252.0 + result * t;
              result = -113.0 / 252.0 + result * t;
              result = -61.0 / 504.0 + result * t;
              result = 1543.0 / 2520.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 11.0 / 280.0;
              result = -13.0 / 84.0 + result * t;
              result = 1.0 / 7.0 + result * t;
              result = 4.0 / 21.0 + result * t;
              result = -3.0 / 7.0 + result * t;
              result = 23.0 / 105.0 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.0 / 120.0;
              result = 1.0 / 24.0 + result * t;
              result = -1.0 / 12.0 + result * t;
              result = 1.0 / 12.0 + result * t;
              result = -1.0 / 24.0 + result * t;
              result = 1.0 / 120.0 + result * t;
              return result;
            }
          }
        }
        return 0.0;
        break;

      default:
        return 0.0;
    }
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override { return bsplineBasis.getDegree(); }

  inline double getIntegral(LT level, IT index) override { return -1.0; }

 protected:
  /// B-spline basis for B-spline evaluation
  BsplineBasis<LT, IT> bsplineBasis;
};

// default type-def (unsigned int for level and index)
typedef NakBsplineBoundaryCombigridBasis<unsigned int, unsigned int>
    SNakBsplineBoundaryCombigridBase;

}  // namespace base
}  // namespace sgpp
