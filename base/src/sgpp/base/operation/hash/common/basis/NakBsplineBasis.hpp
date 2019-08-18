// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NAK_BSPLINE_BASE_HPP
#define NAK_BSPLINE_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <cstddef>
#include <algorithm>

namespace sgpp {
namespace base {

/**
 * Not-a-knot B-spline basis.
 */
template <class LT, class IT>
class NakBsplineBasis: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakBsplineBasis() : bsplineBasis(BsplineBasis<LT, IT>()) {
  }

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineBasis(size_t degree) :
        bsplineBasis(BsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakBsplineBasis() override {
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
          if (i == 0) {
            // l = 1, i = 0
            return 0.5 * t * t - 1.5 * t + 1.0;
          } else if (i == 1) {
            // l = 1, i = 1
            return 1.0 - t * t;
          } else {
            // l = 1, i = 2
            return 0.5 * t * t + 1.5 * t + 1.0;
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

      case 5:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return 1.0 - x;
          } else {
            // l = 0, i = 1
            return x;
          }
        } else if (l == 1) {
          if (i == 0) {
            // l = 1, i = 0
            return 0.5 * t * t - 1.5 * t + 1.0;
          } else if (i == 1) {
            // l = 1, i = 1
            return 1.0 - t * t;
          } else {
            // l = 1, i = 2
            return 0.5 * t * t + 1.5 * t + 1.0;
          }
        } else if ((i > 5) && (i < hInv - 5)) {
          // l >= 4, 5 < i < 2^l - 5
          return bsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if ((l == 2) && (i == 0)) {
            // l = 2, i = 0
            double result = 1.0/24.0;
            result = -5.0/12.0 + result * t;
            result = 35.0/24.0 + result * t;
            result = -25.0/12.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else if ((l == 2) && (i == 1)) {
            // l = 2, i = 1
            double result = -1.0/6.0;
            result = 5.0/6.0 + result * t;
            result = -5.0/6.0 + result * t;
            result = -5.0/6.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else if ((l == 2) && (i == 2)) {
            // l = 2, i = 2
            double result = 1.0/4.0;
            result *= t;
            result = -5.0/4.0 + result * t;
            result *= t;
            result = 1.0 + result * t;
            return result;
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 107.0/30240.0;
              result = -17.0/756.0 + result * t;
              result = 1.0/378.0 + result * t;
              result = 37.0/378.0 + result * t;
              result = 109.0/756.0 + result * t;
              result = 253.0/3780.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -397.0/30240.0;
              result = 185.0/6048.0 + result * t;
              result = 155.0/3024.0 + result * t;
              result = -415.0/3024.0 + result * t;
              result = -1165.0/6048.0 + result * t;
              result = 2965.0/6048.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 233.0/30240.0;
              result = -53.0/1512.0 + result * t;
              result = 8.0/189.0 + result * t;
              result = 13.0/189.0 + result * t;
              result = -97.0/378.0 + result * t;
              result = 433.0/1890.0 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.0/4320.0;
              result = 1.0/288.0 + result * t;
              result = -1.0/48.0 + result * t;
              result = 1.0/16.0 + result * t;
              result = -3.0/32.0 + result * t;
              result = 9.0/160.0 + result * t;
              return result;
            }
          } else if ((l == 3) && (i == 4)) {
            // l = 3, i = 4
            if ((t < -4.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 4.0;
              double result = -1.3888888888888889e-03;
              result = 4.6296296296296294e-03 + result * t;
              result = 9.2592592592592587e-03 + result * t;
              result = 9.2592592592592587e-03 + result * t;
              result = 4.6296296296296294e-03 + result * t;
              result = 9.2592592592592596e-04 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 1.2500000000000001e-02;
              result = -1.6203703703703703e-02 + result * t;
              result = -6.0185185185185182e-02 + result * t;
              result = -3.2407407407407406e-02 + result * t;
              result = 2.4768518518518517e-01 + result * t;
              result = 3.8564814814814813e-01 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -1.2500000000000001e-02;
              result = 4.6296296296296294e-02 + result * t;
              result *= t;
              result = -1.8518518518518517e-01 + result * t;
              result *= t;
              result = 5.3703703703703709e-01 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = 1.3888888888888889e-03;
              result = -1.6203703703703703e-02 + result * t;
              result = 6.0185185185185182e-02 + result * t;
              result = -3.2407407407407406e-02 + result * t;
              result = -2.4768518518518517e-01 + result * t;
              result = 3.8564814814814813e-01 + result * t;
              return result;
            }
          } else if (i == 0) {
            // l >= 3, i = 0
            if ((t < 0.0) || (t > 3.0)) {
              return 0.0;
            } else {
              double result = -3.9682539682539683e-04;
              result = 5.9523809523809521e-03 + result * t;
              result = -3.5714285714285712e-02 + result * t;
              result = 1.0714285714285714e-01 + result * t;
              result = -1.6071428571428573e-01 + result * t;
              result = 9.6428571428571433e-02 + result * t;
              return result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 1.0/504.0;
              result = -1.0/42.0 + result * t;
              result = 2.0/21.0 + result * t;
              result = -2.0/21.0 + result * t;
              result = -5.0/21.0 + result * t;
              result = 47.0/105.0 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.0/840.0;
              result = 1.0/168.0 + result * t;
              result = -1.0/84.0 + result * t;
              result = 1.0/84.0 + result * t;
              result = -1.0/168.0 + result * t;
              result = 1.0/840.0 + result * t;
              return result;
            }
          } else if (i == 2) {
            // l >= 3, i = 2
            if ((t < -2.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 2.0;
              double result = -3.9682539682539680e-03;
              result = 3.5714285714285712e-02 + result * t;
              result = -7.1428571428571425e-02 + result * t;
              result = -1.1904761904761904e-01 + result * t;
              result = 2.5000000000000000e-01 + result * t;
              result = 3.8809523809523810e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 7.1428571428571426e-03;
              result = -2.3809523809523808e-02 + result * t;
              result *= t;
              result = 9.5238095238095233e-02 + result * t;
              result = -1.4285714285714285e-01 + result * t;
              result = 6.6666666666666666e-02 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -2.3809523809523812e-03;
              result = 1.1904761904761904e-02 + result * t;
              result = -2.3809523809523808e-02 + result * t;
              result = 2.3809523809523808e-02 + result * t;
              result = -1.1904761904761904e-02 + result * t;
              result = 2.3809523809523812e-03 + result * t;
              return result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 1.0/252.0;
              result = -1.0/42.0 + result * t;
              result *= t;
              result = 2.0/21.0 + result * t;
              result = 1.0/7.0 + result * t;
              result = 1.0/15.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -23.0/1260.0;
              result = 1.0/28.0 + result * t;
              result = 1.0/14.0 + result * t;
              result = -5.0/42.0 + result * t;
              result = -1.0/4.0 + result * t;
              result = 163.0/420.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 19.0/1260.0;
              result = -1.0/18.0 + result * t;
              result = 2.0/63.0 + result * t;
              result = 8.0/63.0 + result * t;
              result = -2.0/9.0 + result * t;
              result = 34.0/315.0 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.0/252.0;
              result = 5.0/252.0 + result * t;
              result = -5.0/126.0 + result * t;
              result = 5.0/126.0 + result * t;
              result = -5.0/252.0 + result * t;
              result = 1.0/252.0 + result * t;
              return result;
            }
          } else if (i == 4) {
            // l >= 4, i = 4
            if ((t < -4.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 4.0;
              double result = -1.9841269841269840e-03;
              result = 5.9523809523809521e-03 + result * t;
              result = 1.1904761904761904e-02 + result * t;
              result = 1.1904761904761904e-02 + result * t;
              result = 5.9523809523809521e-03 + result * t;
              result = 1.1904761904761906e-03 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 2.5793650793650792e-02;
              result = -2.3809523809523808e-02 + result * t;
              result = -9.5238095238095233e-02 + result * t;
              result = -9.5238095238095233e-02 + result * t;
              result = 2.3809523809523808e-01 + result * t;
              result = 4.4761904761904764e-01 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -4.0873015873015874e-02;
              result = 1.0515873015873016e-01 + result * t;
              result = 6.7460317460317457e-02 + result * t;
              result = -2.6587301587301587e-01 + result * t;
              result = -2.0436507936507936e-01 + result * t;
              result = 4.9722222222222223e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 2.5793650793650792e-02;
              result = -9.9206349206349201e-02 + result * t;
              result = 7.9365079365079361e-02 + result * t;
              result = 1.5873015873015872e-01 + result * t;
              result = -3.1746031746031744e-01 + result * t;
              result = 1.5873015873015872e-01 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -5.9523809523809521e-03;
              result = 2.9761904761904760e-02 + result * t;
              result = -5.9523809523809521e-02 + result * t;
              result = 5.9523809523809521e-02 + result * t;
              result = -2.9761904761904760e-02 + result * t;
              result = 5.9523809523809521e-03 + result * t;
              return result;
            }
          } else {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < -2.0) {
              t += 5.0;
              double result = 1.0/2520.0;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = -11.0/504.0;
              result = 1.0/168.0 + result * t;
              result = 1.0/28.0 + result * t;
              result = 3.0/28.0 + result * t;
              result = 9.0/56.0 + result * t;
              result = 27.0/280.0 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 31.0/504.0;
              result = -13.0/126.0 + result * t;
              result = -10.0/63.0 + result * t;
              result = 2.0/63.0 + result * t;
              result = 25.0/63.0 + result * t;
              result = 121.0/315.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -181.0/2520.0;
              result = 103.0/504.0 + result * t;
              result = 11.0/252.0 + result * t;
              result = -113.0/252.0 + result * t;
              result = -61.0/504.0 + result * t;
              result = 1543.0/2520.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 11.0/280.0;
              result = -13.0/84.0 + result * t;
              result = 1.0/7.0 + result * t;
              result = 4.0/21.0 + result * t;
              result = -3.0/7.0 + result * t;
              result = 23.0/105.0 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.0/120.0;
              result = 1.0/24.0 + result * t;
              result = -1.0/12.0 + result * t;
              result = 1.0/12.0 + result * t;
              result = -1.0/24.0 + result * t;
              result = 1.0/120.0 + result * t;
              return result;
            }
          }
        }

      case 7:
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
        } else if ((i > 7) && (i < hInv - 7)) {
          // l >= 5, 5 < i < 2^l - 5
          return bsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -1.0/6.0;
            result = 5.0/6.0 + result * t;
            result = -5.0/6.0 + result * t;
            result = -5.0/6.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 23.0/4118400.0;
              result = -931.0/6177600.0 + result * t;
              result = 287.0/171600.0 + result * t;
              result = -721.0/77220.0 + result * t;
              result = 49.0/2145.0 + result * t;
              result = 476.0/32175.0 + result * t;
              result = -6608.0/32175.0 + result * t;
              result = 31808.0/96525.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0/4942080.0;
              result = 7.0/1235520.0 + result * t;
              result = -7.0/102960.0 + result * t;
              result = 7.0/15444.0 + result * t;
              result = -7.0/3861.0 + result * t;
              result = 28.0/6435.0 + result * t;
              result = -112.0/19305.0 + result * t;
              result = 64.0/19305.0 + result * t;
              return result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 151.0/12355200.0;
              result = -1589.0/6177600.0 + result * t;
              result = 889.0/514800.0 + result * t;
              result = -119.0/77220.0 + result * t;
              result = -1589.0/77220.0 + result * t;
              result = 161.0/11700.0 + result * t;
              result = 62041.0/386100.0 + result * t;
              result = 68833.0/386100.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -127.0/24710400.0;
              result = 7.0/82368.0 + result * t;
              result = -7.0/20592.0 + result * t;
              result = -7.0/5148.0 + result * t;
              result = 175.0/15444.0 + result * t;
              result = -7.0/8580.0 + result * t;
              result = -2023.0/15444.0 + result * t;
              result = 6299.0/25740.0 + result * t;
              return result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 1.0/86400.0;
              result = -1.0/3600.0 + result * t;
              result = 19.0/7200.0 + result * t;
              result = -17.0/1440.0 + result * t;
              result = 79.0/4320.0 + result * t;
              result = 103.0/2400.0 + result * t;
              result = -4361.0/21600.0 + result * t;
              result = 11023.0/50400.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0/151200.0;
              result = 1.0/21600.0 + result * t;
              result = -1.0/7200.0 + result * t;
              result = 1.0/4320.0 + result * t;
              result = -1.0/4320.0 + result * t;
              result = 1.0/7200.0 + result * t;
              result = -1.0/21600.0 + result * t;
              result = 1.0/151200.0 + result * t;
              return result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 1.0/17280.0;
              result = -1.0/1080.0 + result * t;
              result = 1.0/240.0 + result * t;
              result = 1.0/432.0 + result * t;
              result = -1.0/24.0 + result * t;
              result = -17.0/720.0 + result * t;
              result = 73.0/360.0 + result * t;
              result = 4051.0/15120.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -13.0/60480.0;
              result = 1.0/1440.0 + result * t;
              result = 1.0/720.0 + result * t;
              result = -1.0/144.0 + result * t;
              result = -1.0/216.0 + result * t;
              result = 13.0/240.0 + result * t;
              result = -97.0/1080.0 + result * t;
              result = 83.0/1680.0 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 1.0/6720.0;
              result = -7.0/8640.0 + result * t;
              result = 1.0/960.0 + result * t;
              result = 5.0/1728.0 + result * t;
              result = -7.0/576.0 + result * t;
              result = 53.0/2880.0 + result * t;
              result = -13.0/960.0 + result * t;
              result = 7.0/1728.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0/30240.0;
              result = 1.0/4320.0 + result * t;
              result = -1.0/1440.0 + result * t;
              result = 1.0/864.0 + result * t;
              result = -1.0/864.0 + result * t;
              result = 1.0/1440.0 + result * t;
              result = -1.0/4320.0 + result * t;
              result = 1.0/30240.0 + result * t;
              return result;
            }
          } else if (i == 5) {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 5.0;
              double result = 1.0/28800.0;
              result = -1.0/3600.0 + result * t;
              result = -1.0/7200.0 + result * t;
              result = 1.0/480.0 + result * t;
              result = 29.0/4320.0 + result * t;
              result = 23.0/2400.0 + result * t;
              result = 149.0/21600.0 + result * t;
              result = 103.0/50400.0 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -7.0/9600.0;
              result = 1.0/1440.0 + result * t;
              result = 7.0/1440.0 + result * t;
              result = 1.0/96.0 + result * t;
              result = -23.0/864.0 + result * t;
              result = -19.0/160.0 + result * t;
              result = 217.0/4320.0 + result * t;
              result = 487.0/1120.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 293.0/201600.0;
              result = -127.0/28800.0 + result * t;
              result = -181.0/28800.0 + result * t;
              result = 113.0/5760.0 + result * t;
              result = 899.0/17280.0 + result * t;
              result = -887.0/9600.0 + result * t;
              result = -17461.0/86400.0 + result * t;
              result = 23851.0/67200.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -257.0/201600.0;
              result = 83.0/14400.0 + result * t;
              result = -1.0/450.0 + result * t;
              result = -13.0/480.0 + result * t;
              result = 131.0/4320.0 + result * t;
              result = 199.0/2400.0 + result * t;
              result = -4321.0/21600.0 + result * t;
              result = 6191.0/50400.0 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 47.0/86400.0;
              result = -91.0/28800.0 + result * t;
              result = 161.0/28800.0 + result * t;
              result = 7.0/1920.0 + result * t;
              result = -511.0/17280.0 + result * t;
              result = 469.0/9600.0 + result * t;
              result = -3199.0/86400.0 + result * t;
              result = 323.0/28800.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0/10800.0;
              result = 7.0/10800.0 + result * t;
              result = -7.0/3600.0 + result * t;
              result = 7.0/2160.0 + result * t;
              result = -7.0/2160.0 + result * t;
              result = 7.0/3600.0 + result * t;
              result = -7.0/10800.0 + result * t;
              result = 1.0/10800.0 + result * t;
              return result;
            }
          } else {
            // l >= 4, i = 7
            if ((t < -7.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 1.0/604800.0;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              return result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -47.0/86400.0;
              result = 1.0/21600.0 + result * t;
              result = 1.0/1800.0 + result * t;
              result = 1.0/270.0 + result * t;
              result = 2.0/135.0 + result * t;
              result = 8.0/225.0 + result * t;
              result = 32.0/675.0 + result * t;
              result = 128.0/4725.0 + result * t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 217.0/86400.0;
              result = -13.0/3456.0 + result * t;
              result = -61.0/5760.0 + result * t;
              result = -41.0/3456.0 + result * t;
              result = 59.0/3456.0 + result * t;
              result = 559.0/5760.0 + result * t;
              result = 3059.0/17280.0 + result * t;
              result = 15559.0/120960.0 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -443.0/86400.0;
              result = 199.0/14400.0 + result * t;
              result = 47.0/2400.0 + result * t;
              result = -1.0/30.0 + result * t;
              result = -89.0/720.0 + result * t;
              result = -13.0/400.0 + result * t;
              result = 1141.0/3600.0 + result * t;
              result = 1109.0/2800.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 3499.0/604800.0;
              result = -1907.0/86400.0 + result * t;
              result = -149.0/28800.0 + result * t;
              result = 1597.0/17280.0 + result * t;
              result = 619.0/17280.0 + result * t;
              result = -8867.0/28800.0 + result * t;
              result = -9269.0/86400.0 + result * t;
              result = 333757.0/604800.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -569.0/151200.0;
              result = 199.0/10800.0 + result * t;
              result = -29.0/1800.0 + result * t;
              result = -67.0/1080.0 + result * t;
              result = 31.0/270.0 + result * t;
              result = 167.0/1800.0 + result * t;
              result = -491.0/1350.0 + result * t;
              result = 9203.0/37800.0 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 67.0/50400.0;
              result = -19.0/2400.0 + result * t;
              result = 37.0/2400.0 + result * t;
              result = 1.0/480.0 + result * t;
              result = -83.0/1440.0 + result * t;
              result = 81.0/800.0 + result * t;
              result = -563.0/7200.0 + result * t;
              result = 401.0/16800.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0/5040.0;
              result = 1.0/720.0 + result * t;
              result = -1.0/240.0 + result * t;
              result = 1.0/144.0 + result * t;
              result = -1.0/144.0 + result * t;
              result = 1.0/240.0 + result * t;
              result = -1.0/720.0 + result * t;
              result = 1.0/5040.0 + result * t;
              return result;
            }
          }
        }

      default:
        return 0.0;
    }
  }

  inline double getIntegral(LT level, IT index) override { return -1.0; }

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
typedef NakBsplineBasis<unsigned int, unsigned int> SNakBsplineBase;

}  // namespace base
}  // namespace sgpp

#endif /* NAK_BSPLINE_BASE_HPP */
