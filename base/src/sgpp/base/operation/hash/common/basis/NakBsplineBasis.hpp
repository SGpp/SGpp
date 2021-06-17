// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/tools/Distribution.hpp>
#include <sgpp/base/tools/DistributionUniform.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>

namespace sgpp {
namespace base {

/**
 * Not-a-knot B-spline basis.
 */
template <class LT, class IT>
class NakBsplineBasis : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakBsplineBasis() : bsplineBasis(BsplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineBasis(size_t degree) : bsplineBasis(BsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakBsplineBasis() override {}

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
            // return 1 - 3 * x + 2 * x * x;
            return 0.5 * t * t - 1.5 * t + 1.0;
          } else if (i == 1) {
            // l = 1, i = 1
            // return 4 * x - 4 * x * x;
            return 1.0 - t * t;
          } else {
            // l = 1, i = 2
            // return 2 * x * x - x;
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
            // return 1 - 3 * x + 2 * x * x;
            return 0.5 * t * t - 1.5 * t + 1.0;
          } else if (i == 1) {
            // l = 1, i = 1
            // return 4 * x - 4 * x * x;
            return 1.0 - t * t;
          } else {
            // l = 1, i = 2
            // return 2 * x * x - x;
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
            double result = 1.0 / 24.0;
            result = -5.0 / 12.0 + result * t;
            result = 35.0 / 24.0 + result * t;
            result = -25.0 / 12.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else if ((l == 2) && (i == 1)) {
            // l = 2, i = 1
            double result = -1.0 / 6.0;
            result = 5.0 / 6.0 + result * t;
            result = -5.0 / 6.0 + result * t;
            result = -5.0 / 6.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else if ((l == 2) && (i == 2)) {
            // l = 2, i = 2
            double result = 1.0 / 4.0;
            result *= t;
            result = -5.0 / 4.0 + result * t;
            result *= t;
            result = 1.0 + result * t;
            return result;
          } else if ((l == 3) && (i == 3)) {
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
            double result = -1.0 / 6.0;
            result = 5.0 / 6.0 + result * t;
            result = -5.0 / 6.0 + result * t;
            result = -5.0 / 6.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 23.0 / 4118400.0;
              result = -931.0 / 6177600.0 + result * t;
              result = 287.0 / 171600.0 + result * t;
              result = -721.0 / 77220.0 + result * t;
              result = 49.0 / 2145.0 + result * t;
              result = 476.0 / 32175.0 + result * t;
              result = -6608.0 / 32175.0 + result * t;
              result = 31808.0 / 96525.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0 / 4942080.0;
              result = 7.0 / 1235520.0 + result * t;
              result = -7.0 / 102960.0 + result * t;
              result = 7.0 / 15444.0 + result * t;
              result = -7.0 / 3861.0 + result * t;
              result = 28.0 / 6435.0 + result * t;
              result = -112.0 / 19305.0 + result * t;
              result = 64.0 / 19305.0 + result * t;
              return result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 151.0 / 12355200.0;
              result = -1589.0 / 6177600.0 + result * t;
              result = 889.0 / 514800.0 + result * t;
              result = -119.0 / 77220.0 + result * t;
              result = -1589.0 / 77220.0 + result * t;
              result = 161.0 / 11700.0 + result * t;
              result = 62041.0 / 386100.0 + result * t;
              result = 68833.0 / 386100.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -127.0 / 24710400.0;
              result = 7.0 / 82368.0 + result * t;
              result = -7.0 / 20592.0 + result * t;
              result = -7.0 / 5148.0 + result * t;
              result = 175.0 / 15444.0 + result * t;
              result = -7.0 / 8580.0 + result * t;
              result = -2023.0 / 15444.0 + result * t;
              result = 6299.0 / 25740.0 + result * t;
              return result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 1.0 / 86400.0;
              result = -1.0 / 3600.0 + result * t;
              result = 19.0 / 7200.0 + result * t;
              result = -17.0 / 1440.0 + result * t;
              result = 79.0 / 4320.0 + result * t;
              result = 103.0 / 2400.0 + result * t;
              result = -4361.0 / 21600.0 + result * t;
              result = 11023.0 / 50400.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0 / 151200.0;
              result = 1.0 / 21600.0 + result * t;
              result = -1.0 / 7200.0 + result * t;
              result = 1.0 / 4320.0 + result * t;
              result = -1.0 / 4320.0 + result * t;
              result = 1.0 / 7200.0 + result * t;
              result = -1.0 / 21600.0 + result * t;
              result = 1.0 / 151200.0 + result * t;
              return result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 1.0 / 17280.0;
              result = -1.0 / 1080.0 + result * t;
              result = 1.0 / 240.0 + result * t;
              result = 1.0 / 432.0 + result * t;
              result = -1.0 / 24.0 + result * t;
              result = -17.0 / 720.0 + result * t;
              result = 73.0 / 360.0 + result * t;
              result = 4051.0 / 15120.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -13.0 / 60480.0;
              result = 1.0 / 1440.0 + result * t;
              result = 1.0 / 720.0 + result * t;
              result = -1.0 / 144.0 + result * t;
              result = -1.0 / 216.0 + result * t;
              result = 13.0 / 240.0 + result * t;
              result = -97.0 / 1080.0 + result * t;
              result = 83.0 / 1680.0 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 1.0 / 6720.0;
              result = -7.0 / 8640.0 + result * t;
              result = 1.0 / 960.0 + result * t;
              result = 5.0 / 1728.0 + result * t;
              result = -7.0 / 576.0 + result * t;
              result = 53.0 / 2880.0 + result * t;
              result = -13.0 / 960.0 + result * t;
              result = 7.0 / 1728.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0 / 30240.0;
              result = 1.0 / 4320.0 + result * t;
              result = -1.0 / 1440.0 + result * t;
              result = 1.0 / 864.0 + result * t;
              result = -1.0 / 864.0 + result * t;
              result = 1.0 / 1440.0 + result * t;
              result = -1.0 / 4320.0 + result * t;
              result = 1.0 / 30240.0 + result * t;
              return result;
            }
          } else if (i == 5) {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 5.0;
              double result = 1.0 / 28800.0;
              result = -1.0 / 3600.0 + result * t;
              result = -1.0 / 7200.0 + result * t;
              result = 1.0 / 480.0 + result * t;
              result = 29.0 / 4320.0 + result * t;
              result = 23.0 / 2400.0 + result * t;
              result = 149.0 / 21600.0 + result * t;
              result = 103.0 / 50400.0 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -7.0 / 9600.0;
              result = 1.0 / 1440.0 + result * t;
              result = 7.0 / 1440.0 + result * t;
              result = 1.0 / 96.0 + result * t;
              result = -23.0 / 864.0 + result * t;
              result = -19.0 / 160.0 + result * t;
              result = 217.0 / 4320.0 + result * t;
              result = 487.0 / 1120.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 293.0 / 201600.0;
              result = -127.0 / 28800.0 + result * t;
              result = -181.0 / 28800.0 + result * t;
              result = 113.0 / 5760.0 + result * t;
              result = 899.0 / 17280.0 + result * t;
              result = -887.0 / 9600.0 + result * t;
              result = -17461.0 / 86400.0 + result * t;
              result = 23851.0 / 67200.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -257.0 / 201600.0;
              result = 83.0 / 14400.0 + result * t;
              result = -1.0 / 450.0 + result * t;
              result = -13.0 / 480.0 + result * t;
              result = 131.0 / 4320.0 + result * t;
              result = 199.0 / 2400.0 + result * t;
              result = -4321.0 / 21600.0 + result * t;
              result = 6191.0 / 50400.0 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 47.0 / 86400.0;
              result = -91.0 / 28800.0 + result * t;
              result = 161.0 / 28800.0 + result * t;
              result = 7.0 / 1920.0 + result * t;
              result = -511.0 / 17280.0 + result * t;
              result = 469.0 / 9600.0 + result * t;
              result = -3199.0 / 86400.0 + result * t;
              result = 323.0 / 28800.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0 / 10800.0;
              result = 7.0 / 10800.0 + result * t;
              result = -7.0 / 3600.0 + result * t;
              result = 7.0 / 2160.0 + result * t;
              result = -7.0 / 2160.0 + result * t;
              result = 7.0 / 3600.0 + result * t;
              result = -7.0 / 10800.0 + result * t;
              result = 1.0 / 10800.0 + result * t;
              return result;
            }
          } else {
            // l >= 4, i = 7
            if ((t < -7.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 1.0 / 604800.0;
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
              double result = -47.0 / 86400.0;
              result = 1.0 / 21600.0 + result * t;
              result = 1.0 / 1800.0 + result * t;
              result = 1.0 / 270.0 + result * t;
              result = 2.0 / 135.0 + result * t;
              result = 8.0 / 225.0 + result * t;
              result = 32.0 / 675.0 + result * t;
              result = 128.0 / 4725.0 + result * t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 217.0 / 86400.0;
              result = -13.0 / 3456.0 + result * t;
              result = -61.0 / 5760.0 + result * t;
              result = -41.0 / 3456.0 + result * t;
              result = 59.0 / 3456.0 + result * t;
              result = 559.0 / 5760.0 + result * t;
              result = 3059.0 / 17280.0 + result * t;
              result = 15559.0 / 120960.0 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -443.0 / 86400.0;
              result = 199.0 / 14400.0 + result * t;
              result = 47.0 / 2400.0 + result * t;
              result = -1.0 / 30.0 + result * t;
              result = -89.0 / 720.0 + result * t;
              result = -13.0 / 400.0 + result * t;
              result = 1141.0 / 3600.0 + result * t;
              result = 1109.0 / 2800.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 3499.0 / 604800.0;
              result = -1907.0 / 86400.0 + result * t;
              result = -149.0 / 28800.0 + result * t;
              result = 1597.0 / 17280.0 + result * t;
              result = 619.0 / 17280.0 + result * t;
              result = -8867.0 / 28800.0 + result * t;
              result = -9269.0 / 86400.0 + result * t;
              result = 333757.0 / 604800.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -569.0 / 151200.0;
              result = 199.0 / 10800.0 + result * t;
              result = -29.0 / 1800.0 + result * t;
              result = -67.0 / 1080.0 + result * t;
              result = 31.0 / 270.0 + result * t;
              result = 167.0 / 1800.0 + result * t;
              result = -491.0 / 1350.0 + result * t;
              result = 9203.0 / 37800.0 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 67.0 / 50400.0;
              result = -19.0 / 2400.0 + result * t;
              result = 37.0 / 2400.0 + result * t;
              result = 1.0 / 480.0 + result * t;
              result = -83.0 / 1440.0 + result * t;
              result = 81.0 / 800.0 + result * t;
              result = -563.0 / 7200.0 + result * t;
              result = 401.0 / 16800.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0 / 5040.0;
              result = 1.0 / 720.0 + result * t;
              result = -1.0 / 240.0 + result * t;
              result = 1.0 / 144.0 + result * t;
              result = -1.0 / 144.0 + result * t;
              result = 1.0 / 240.0 + result * t;
              result = -1.0 / 720.0 + result * t;
              result = 1.0 / 5040.0 + result * t;
              return result;
            }
          }
        }

      default:
        return 0.0;
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of derivative of B-spline basis function
   */
  inline double evalDx(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    double innerDeriv = hInvDbl;
    double t = x * hInvDbl - static_cast<double>(i);

    switch (getDegree()) {
      case 1:
        if ((t < -2.0) || (t > 2.0)) {
          return 0.0;
        } else if (t < 0.0) {
          return innerDeriv;
        } else {
          return -innerDeriv;
        }

      case 3:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return -innerDeriv;
          } else {
            // l = 0, i = 1
            return innerDeriv;
          }
        } else if (l == 1) {
          if (i == 0) {
            // l = 1, i = 0
            return innerDeriv * (t - 1.5);
          } else if (i == 1) {
            // l = 1, i = 1
            return innerDeriv * (-2.0 * t);
          } else {
            // l = 1, i = 2
            return innerDeriv * (t + 1.5);
          }
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 4, 3 < i < 2^l - 3
          return bsplineBasis.evalDx(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }

          if (i == 0) {
            // l >= 2, i = 0
            if ((t < 0.0) || (t > 2.0)) {
              return 0.0;
            } else {
              double result = -1.2500000000000000e-01;
              result = 5.0000000000000000e-01 + result * t;
              result = -5.0000000000000000e-01 + result * t;
              return innerDeriv * result;
            }
          } else if ((l == 2) && (i == 1)) {
            // l = 2, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 2.9999999999999999e-01;
              result = -9.0000000000000002e-01 + result * t;
              result = 2.9999999999999999e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -7.4999999999999997e-02;
              result = 2.9999999999999999e-01 + result * t;
              result = -2.9999999999999999e-01 + result * t;
              return innerDeriv * result;
            }
          } else if (l == 2) {
            // l = 2, i = 2
            if ((t < -2.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 2.0;
              double result = -2.5000000000000000e-01;
              result = 4.0000000000000002e-01 + result * t;
              result = 2.0000000000000001e-01 + result * t;
              return innerDeriv * result;
            } else {
              double result = 2.5000000000000000e-01;
              result = -5.9999999999999998e-01 + result * t;
              result *= t;
              return innerDeriv * result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 3.7500000000000000e-01;
              result = -1.0 + result * t;
              result = 2.5000000000000000e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -2.5000000000000000e-01;
              result = 5.0000000000000000e-01 + result * t;
              result = -2.5000000000000000e-01 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 2) {
            // l >= 3, i = 2
            if ((t < -2.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 2.0;
              double result = -3.7500000000000000e-01;
              result = 5.0000000000000000e-01 + result * t;
              result = 2.5000000000000000e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 8.7500000000000000e-01;
              result = -1.0 + result * t;
              result = -2.5000000000000000e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -3.7500000000000000e-01;
              result = 7.5000000000000000e-01 + result * t;
              result = -3.7500000000000000e-01 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 3, i = 3
            if ((t < -3.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 3.0;
              double result = 1.2500000000000000e-01;
              result *= t;
              result *= t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -1.1250000000000000e+00;
              result = 5.0000000000000000e-01 + result * t;
              result = 5.0000000000000000e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 1.3750000000000000e+00;
              result = -1.7500000000000000e+00 + result * t;
              result = -1.2500000000000000e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -5.0000000000000000e-01;
              result = 1.0 + result * t;
              result = -5.0000000000000000e-01 + result * t;
              return innerDeriv * result;
            }
          }
        }

      case 5:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return -innerDeriv;
          } else {
            // l = 0, i = 1
            return innerDeriv;
          }
        } else if (l == 1) {
          if (i == 0) {
            // l = 1, i = 0
            return innerDeriv * (t - 1.5);
          } else if (i == 1) {
            // l = 1, i = 1
            return innerDeriv * (-2.0 * t);
          } else {
            // l = 1, i = 2
            return innerDeriv * (t + 1.5);
          }
        } else if ((i > 5) && (i < hInv - 5)) {
          // l >= 4, 5 < i < 2^l - 5
          return bsplineBasis.evalDx(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }

          if ((l == 2) && (i == 0)) {
            // l = 2, i = 0
            double result = 1.0 / 6.0;
            result = -5.0 / 4.0 + result * t;
            result = 35.0 / 12.0 + result * t;
            result = -25.0 / 12.0 + result * t;
            return innerDeriv * result;
          } else if ((l == 2) && (i == 1)) {
            // l = 2, i = 1
            double result = -2.0 / 3.0;
            result = 5.0 / 2.0 + result * t;
            result = -5.0 / 3.0 + result * t;
            result = -5.0 / 6.0 + result * t;
            return innerDeriv * result;
          } else if ((l == 2) && (i == 2)) {
            // l = 2, i = 2
            double result = 1.0;
            result *= t;
            result = -5.0 / 2.0 + result * t;
            result *= t;
            return innerDeriv * result;
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 1.7691798941798943e-02;
              result = -8.9947089947089942e-02 + result * t;
              result = 7.9365079365079361e-03 + result * t;
              result = 1.9576719576719576e-01 + result * t;
              result = 1.4417989417989419e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -6.5641534391534390e-02;
              result = 1.2235449735449735e-01 + result * t;
              result = 1.5376984126984128e-01 + result * t;
              result = -2.7447089947089948e-01 + result * t;
              result = -1.9262566137566137e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 3.8525132275132275e-02;
              result = -1.4021164021164020e-01 + result * t;
              result = 1.2698412698412698e-01 + result * t;
              result = 1.3756613756613756e-01 + result * t;
              result = -2.5661375661375663e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -1.1574074074074073e-03;
              result = 1.3888888888888888e-02 + result * t;
              result = -6.2500000000000000e-02 + result * t;
              result = 1.2500000000000000e-01 + result * t;
              result = -9.3750000000000000e-02 + result * t;
              return innerDeriv * result;
            }
          } else if ((l == 3) && (i == 4)) {
            // l = 3, i = 4
            if ((t < -4.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 4.0;
              double result = -6.9444444444444441e-03;
              result = 1.8518518518518517e-02 + result * t;
              result = 2.7777777777777776e-02 + result * t;
              result = 1.8518518518518517e-02 + result * t;
              result = 4.6296296296296294e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 6.2500000000000000e-02;
              result = -6.4814814814814811e-02 + result * t;
              result = -1.8055555555555555e-01 + result * t;
              result = -6.4814814814814811e-02 + result * t;
              result = 2.4768518518518517e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -6.2500000000000000e-02;
              result = 1.8518518518518517e-01 + result * t;
              result *= t;
              result = -3.7037037037037035e-01 + result * t;
              result *= t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = 6.9444444444444441e-03;
              result = -6.4814814814814811e-02 + result * t;
              result = 1.8055555555555555e-01 + result * t;
              result = -6.4814814814814811e-02 + result * t;
              result = -2.4768518518518517e-01 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 0) {
            // l >= 3, i = 0
            if ((t < 0.0) || (t > 3.0)) {
              return 0.0;
            } else {
              double result = -1.9841269841269840e-03;
              result = 2.3809523809523808e-02 + result * t;
              result = -1.0714285714285714e-01 + result * t;
              result = 2.1428571428571427e-01 + result * t;
              result = -1.6071428571428573e-01 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 9.9206349206349201e-03;
              result = -9.5238095238095233e-02 + result * t;
              result = 2.8571428571428570e-01 + result * t;
              result = -1.9047619047619047e-01 + result * t;
              result = -2.3809523809523808e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -5.9523809523809521e-03;
              result = 2.3809523809523808e-02 + result * t;
              result = -3.5714285714285712e-02 + result * t;
              result = 2.3809523809523808e-02 + result * t;
              result = -5.9523809523809521e-03 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 2) {
            // l >= 3, i = 2
            if ((t < -2.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 2.0;
              double result = -1.9841269841269840e-02;
              result = 1.4285714285714285e-01 + result * t;
              result = -2.1428571428571427e-01 + result * t;
              result = -2.3809523809523808e-01 + result * t;
              result = 2.5000000000000000e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 3.5714285714285712e-02;
              result = -9.5238095238095233e-02 + result * t;
              result *= t;
              result = 1.9047619047619047e-01 + result * t;
              result = -1.4285714285714285e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -1.1904761904761904e-02;
              result = 4.7619047619047616e-02 + result * t;
              result = -7.1428571428571425e-02 + result * t;
              result = 4.7619047619047616e-02 + result * t;
              result = -1.1904761904761904e-02 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 1.9841269841269840e-02;
              result = -9.5238095238095233e-02 + result * t;
              result *= t;
              result = 1.9047619047619047e-01 + result * t;
              result = 1.4285714285714285e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -9.1269841269841265e-02;
              result = 1.4285714285714285e-01 + result * t;
              result = 2.1428571428571427e-01 + result * t;
              result = -2.3809523809523808e-01 + result * t;
              result = -2.5000000000000000e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 7.5396825396825393e-02;
              result = -2.2222222222222221e-01 + result * t;
              result = 9.5238095238095233e-02 + result * t;
              result = 2.5396825396825395e-01 + result * t;
              result = -2.2222222222222221e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -1.9841269841269840e-02;
              result = 7.9365079365079361e-02 + result * t;
              result = -1.1904761904761904e-01 + result * t;
              result = 7.9365079365079361e-02 + result * t;
              result = -1.9841269841269840e-02 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 4) {
            // l >= 4, i = 4
            if ((t < -4.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 4.0;
              double result = -9.9206349206349201e-03;
              result = 2.3809523809523808e-02 + result * t;
              result = 3.5714285714285712e-02 + result * t;
              result = 2.3809523809523808e-02 + result * t;
              result = 5.9523809523809521e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 1.2896825396825398e-01;
              result = -9.5238095238095233e-02 + result * t;
              result = -2.8571428571428570e-01 + result * t;
              result = -1.9047619047619047e-01 + result * t;
              result = 2.3809523809523808e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -2.0436507936507936e-01;
              result = 4.2063492063492064e-01 + result * t;
              result = 2.0238095238095238e-01 + result * t;
              result = -5.3174603174603174e-01 + result * t;
              result = -2.0436507936507936e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 1.2896825396825398e-01;
              result = -3.9682539682539680e-01 + result * t;
              result = 2.3809523809523808e-01 + result * t;
              result = 3.1746031746031744e-01 + result * t;
              result = -3.1746031746031744e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -2.9761904761904760e-02;
              result = 1.1904761904761904e-01 + result * t;
              result = -1.7857142857142858e-01 + result * t;
              result = 1.1904761904761904e-01 + result * t;
              result = -2.9761904761904760e-02 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < -2.0) {
              t += 5.0;
              double result = 1.9841269841269840e-03;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              return innerDeriv * result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = -1.0912698412698413e-01;
              result = 2.3809523809523808e-02 + result * t;
              result = 1.0714285714285714e-01 + result * t;
              result = 2.1428571428571427e-01 + result * t;
              result = 1.6071428571428573e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 3.0753968253968256e-01;
              result = -4.1269841269841268e-01 + result * t;
              result = -4.7619047619047616e-01 + result * t;
              result = 6.3492063492063489e-02 + result * t;
              result = 3.9682539682539680e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -3.5912698412698413e-01;
              result = 8.1746031746031744e-01 + result * t;
              result = 1.3095238095238096e-01 + result * t;
              result = -8.9682539682539686e-01 + result * t;
              result = -1.2103174603174603e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 1.9642857142857142e-01;
              result = -6.1904761904761907e-01 + result * t;
              result = 4.2857142857142855e-01 + result * t;
              result = 3.8095238095238093e-01 + result * t;
              result = -4.2857142857142855e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -4.1666666666666664e-02;
              result = 1.6666666666666666e-01 + result * t;
              result = -2.5000000000000000e-01 + result * t;
              result = 1.6666666666666666e-01 + result * t;
              result = -4.1666666666666664e-02 + result * t;
              return innerDeriv * result;
            }
          }
        }

      case 7:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return -innerDeriv;
          } else {
            // l = 0, i = 1
            return innerDeriv;
          }
        } else if (l == 1) {
          // l = 1, i = 1
          return innerDeriv * (-2.0 * t);
        } else if ((i > 7) && (i < hInv - 7)) {
          // l >= 5, 5 < i < 2^l - 5
          return bsplineBasis.evalDx(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -2.0 / 3.0;
            result = 5.0 / 2.0 + result * t;
            result = -5.0 / 3.0 + result * t;
            result = -5.0 / 6.0 + result * t;
            return innerDeriv * result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 3.9092851592851593e-05;
              result = -9.0423465423465426e-04 + result * t;
              result = 8.3624708624708624e-03 + result * t;
              result = -3.7347837347837351e-02 + result * t;
              result = 6.8531468531468534e-02 + result * t;
              result = 2.9588189588189588e-02 + result * t;
              result = -2.0537684537684539e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -1.4164076664076663e-06;
              result = 3.3993783993783994e-05 + result * t;
              result = -3.3993783993783994e-04 + result * t;
              result = 1.8130018130018131e-03 + result * t;
              result = -5.4390054390054390e-03 + result * t;
              result = 8.7024087024087024e-03 + result * t;
              result = -5.8016058016058013e-03 + result * t;
              return innerDeriv * result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 8.5551023051023053e-05;
              result = -1.5433177933177934e-03 + result * t;
              result = 8.6344211344211337e-03 + result * t;
              result = -6.1642061642061645e-03 + result * t;
              result = -6.1732711732711734e-02 + result * t;
              result = 2.7521367521367520e-02 + result * t;
              result = 1.6068635068635068e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -3.5976754726754728e-05;
              result = 5.0990675990675988e-04 + result * t;
              result = -1.6996891996891997e-03 + result * t;
              result = -5.4390054390054390e-03 + result * t;
              result = 3.3993783993783992e-02 + result * t;
              result = -1.6317016317016317e-03 + result * t;
              result = -1.3098938098938098e-01 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 8.1018518518518516e-05;
              result = -1.6666666666666668e-03 + result * t;
              result = 1.3194444444444444e-02 + result * t;
              result = -4.7222222222222221e-02 + result * t;
              result = 5.4861111111111110e-02 + result * t;
              result = 8.5833333333333331e-02 + result * t;
              result = -2.0189814814814816e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -4.6296296296296294e-05;
              result = 2.7777777777777778e-04 + result * t;
              result = -6.9444444444444447e-04 + result * t;
              result = 9.2592592592592596e-04 + result * t;
              result = -6.9444444444444447e-04 + result * t;
              result = 2.7777777777777778e-04 + result * t;
              result = -4.6296296296296294e-05 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 4.0509259259259258e-04;
              result = -5.5555555555555558e-03 + result * t;
              result = 2.0833333333333332e-02 + result * t;
              result = 9.2592592592592587e-03 + result * t;
              result = -1.2500000000000000e-01 + result * t;
              result = -4.7222222222222221e-02 + result * t;
              result = 2.0277777777777778e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -1.5046296296296296e-03;
              result = 4.1666666666666666e-03 + result * t;
              result = 6.9444444444444441e-03 + result * t;
              result = -2.7777777777777776e-02 + result * t;
              result = -1.3888888888888888e-02 + result * t;
              result = 1.0833333333333334e-01 + result * t;
              result = -8.9814814814814820e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 1.0416666666666667e-03;
              result = -4.8611111111111112e-03 + result * t;
              result = 5.2083333333333330e-03 + result * t;
              result = 1.1574074074074073e-02 + result * t;
              result = -3.6458333333333336e-02 + result * t;
              result = 3.6805555555555557e-02 + result * t;
              result = -1.3541666666666667e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -2.3148148148148149e-04;
              result = 1.3888888888888889e-03 + result * t;
              result = -3.4722222222222220e-03 + result * t;
              result = 4.6296296296296294e-03 + result * t;
              result = -3.4722222222222220e-03 + result * t;
              result = 1.3888888888888889e-03 + result * t;
              result = -2.3148148148148149e-04 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 5) {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 5.0;
              double result = 2.4305555555555555e-04;
              result = -1.6666666666666668e-03 + result * t;
              result = -6.9444444444444447e-04 + result * t;
              result = 8.3333333333333332e-03 + result * t;
              result = 2.0138888888888890e-02 + result * t;
              result = 1.9166666666666665e-02 + result * t;
              result = 6.8981481481481480e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -5.1041666666666666e-03;
              result = 4.1666666666666666e-03 + result * t;
              result = 2.4305555555555556e-02 + result * t;
              result = 4.1666666666666664e-02 + result * t;
              result = -7.9861111111111105e-02 + result * t;
              result = -2.3749999999999999e-01 + result * t;
              result = 5.0231481481481481e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 1.0173611111111111e-02;
              result = -2.6458333333333334e-02 + result * t;
              result = -3.1423611111111110e-02 + result * t;
              result = 7.8472222222222221e-02 + result * t;
              result = 1.5607638888888889e-01 + result * t;
              result = -1.8479166666666666e-01 + result * t;
              result = -2.0209490740740740e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -8.9236111111111113e-03;
              result = 3.4583333333333334e-02 + result * t;
              result = -1.1111111111111112e-02 + result * t;
              result = -1.0833333333333334e-01 + result * t;
              result = 9.0972222222222218e-02 + result * t;
              result = 1.6583333333333333e-01 + result * t;
              result = -2.0004629629629631e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 3.8078703703703703e-03;
              result = -1.8958333333333334e-02 + result * t;
              result = 2.7951388888888890e-02 + result * t;
              result = 1.4583333333333334e-02 + result * t;
              result = -8.8715277777777782e-02 + result * t;
              result = 9.7708333333333328e-02 + result * t;
              result = -3.7025462962962961e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -6.4814814814814813e-04;
              result = 3.8888888888888888e-03 + result * t;
              result = -9.7222222222222224e-03 + result * t;
              result = 1.2962962962962963e-02 + result * t;
              result = -9.7222222222222224e-03 + result * t;
              result = 3.8888888888888888e-03 + result * t;
              result = -6.4814814814814813e-04 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 7
            if ((t < -7.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 1.1574074074074073e-05;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              return innerDeriv * result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -3.8078703703703703e-03;
              result = 2.7777777777777778e-04 + result * t;
              result = 2.7777777777777779e-03 + result * t;
              result = 1.4814814814814815e-02 + result * t;
              result = 4.4444444444444446e-02 + result * t;
              result = 7.1111111111111111e-02 + result * t;
              result = 4.7407407407407405e-02 + result * t;
              return innerDeriv * result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 1.7581018518518520e-02;
              result = -2.2569444444444444e-02 + result * t;
              result = -5.2951388888888888e-02 + result * t;
              result = -4.7453703703703706e-02 + result * t;
              result = 5.1215277777777776e-02 + result * t;
              result = 1.9409722222222223e-01 + result * t;
              result = 1.7702546296296295e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -3.5891203703703703e-02;
              result = 8.2916666666666666e-02 + result * t;
              result = 9.7916666666666666e-02 + result * t;
              result = -1.3333333333333333e-01 + result * t;
              result = -3.7083333333333335e-01 + result * t;
              result = -6.5000000000000002e-02 + result * t;
              result = 3.1694444444444442e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 4.0497685185185185e-02;
              result = -1.3243055555555555e-01 + result * t;
              result = -2.5868055555555554e-02 + result * t;
              result = 3.6967592592592591e-01 + result * t;
              result = 1.0746527777777778e-01 + result * t;
              result = -6.1576388888888889e-01 + result * t;
              result = -1.0728009259259259e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -2.6342592592592591e-02;
              result = 1.1055555555555556e-01 + result * t;
              result = -8.0555555555555561e-02 + result * t;
              result = -2.4814814814814815e-01 + result * t;
              result = 3.4444444444444444e-01 + result * t;
              result = 1.8555555555555556e-01 + result * t;
              result = -3.6370370370370370e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 9.3055555555555548e-03;
              result = -4.7500000000000001e-02 + result * t;
              result = 7.7083333333333337e-02 + result * t;
              result = 8.3333333333333332e-03 + result * t;
              result = -1.7291666666666666e-01 + result * t;
              result = 2.0250000000000001e-01 + result * t;
              result = -7.8194444444444441e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -1.3888888888888889e-03;
              result = 8.3333333333333332e-03 + result * t;
              result = -2.0833333333333332e-02 + result * t;
              result = 2.7777777777777776e-02 + result * t;
              result = -2.0833333333333332e-02 + result * t;
              result = 8.3333333333333332e-03 + result * t;
              result = -1.3888888888888889e-03 + result * t;
              return innerDeriv * result;
            }
          }
        }

      default:
        return 0.0;
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @return      integral of basis function
   */
  inline double getIntegral(LT l, IT i) override {
    size_t quadOrder = getDegree() + 1;
    auto pdf_uniform = std::make_shared<sgpp::base::DistributionUniform>(0, 1);
    base::DataVector temp_quadCoordinates, temp_quadWeights;
    base::GaussLegendreQuadRule1D gauss;
    gauss.getLevelPointsAndWeightsNormalized(quadOrder, temp_quadCoordinates, temp_quadWeights);
    auto quadCoordinates = std::make_shared<sgpp::base::DataVector>(temp_quadCoordinates);
    auto quadWeights = std::make_shared<sgpp::base::DataVector>(temp_quadWeights);
    return getMean(l, i, pdf_uniform, quadCoordinates, quadWeights);
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param pdf   probability density function
   * @param quadCoordinates coordinates of the quadrature rule to be used
   * @param quadWeights weights of the quadrature rule to be used
   * @return      integral of basis function
   */
  inline double getMean(LT l, IT i, std::shared_ptr<sgpp::base::Distribution> pdf,
                        std::shared_ptr<sgpp::base::DataVector> quadCoordinates,
                        std::shared_ptr<sgpp::base::DataVector> quadWeights) {
    size_t degree = getDegree();
    if ((degree != 1) && (degree != 3) && (degree != 5)) {
      throw std::runtime_error(
          "NakBsplineBasis: only B spline degrees 1, 3 and 5 are "
          "supported.");
    }

    const size_t pp1h = (degree + 1) >> 1;  //  =|_(p+1)/2_|
    const size_t hInv = 1 << l;             // = 2^lid
    const double hik = 1.0 / static_cast<double>(hInv);
    double offset = (i - static_cast<double>(pp1h)) * hik;

    if (degree == 3) {
      if (i == 3) offset -= hik;
    } else if (degree == 5) {
      if (i == 5) offset -= 2 * hik;
    }

    // start and stop identify the segments on which the spline is nonzero
    size_t start = 0, stop = 0;
    start = ((i > pp1h) ? 0 : (pp1h - i));
    stop = std::min(degree, hInv + pp1h - i - 1);
    // nak special cases
    if (degree == 3) {
      if ((i == 3) || (i == hInv - 3)) stop += 1;
    } else if (degree == 5) {
      if ((i == 5) || (i == hInv - 5)) stop += 2;
    }
    if (l == 2) {
      start = 1;
      stop = 4;
      offset = -0.25;
    }
    if ((degree == 5) && (l == 3)) {
      start = 1;
      stop = 8;
      offset = -0.125;
    }

    double temp_res = basisMean(l, i, start, stop, offset, hik, quadCoordinates, quadWeights, pdf);
    double mean = temp_res * hik;

    return mean;
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override { return bsplineBasis.getDegree(); }

 protected:
  /// B-spline basis for B-spline evaluation
  BsplineBasis<LT, IT> bsplineBasis;

 private:
  /**
   * integrate one basis function
   *
   * @param l				level
   * @param i				index
   * @param start			index for the supports first segment (usually 0)
   * @param stop			index for the supports last segment (usually degree)
   * @param offset			left point of the support
   * @param hik       grid width
   * @param quadCoordinates	the quadrature points
   * @param quadWeights		the quadrature weights
   * @param pdf probability density function
   *
   * @return integral
   */
  double basisMean(LT l, IT i, size_t start, size_t stop, double offset, double hik,
                   std::shared_ptr<base::DataVector> quadCoordinates,
                   std::shared_ptr<base::DataVector> quadWeights,
                   std::shared_ptr<sgpp::base::Distribution> pdf) {
    sgpp::base::DataVector bounds = pdf->getBounds();
    double left = bounds[0];
    double right = bounds[1];

    double temp_res = 0.0;
    // loop over the segments the B-spline is defined on
    for (size_t n = start; n <= stop; n++) {
      // loop over quadrature points
      for (size_t c = 0; c < quadCoordinates->getSize(); c++) {
        // transform  the quadrature points to the segment on which the Bspline is
        // evaluated and the support of the pdf
        const double x = offset + hik * (quadCoordinates->get(c) + static_cast<double>(n));
        double scaledX = left + (right - left) * x;
        temp_res += quadWeights->get(c) * this->eval(l, i, x) * pdf->eval(scaledX);
      }
    }
    return temp_res * (right - left);
  }
};

// default type-def (unsigned int for level and index)
typedef NakBsplineBasis<unsigned int, unsigned int> SNakBsplineBase;

}  // namespace base
}  // namespace sgpp
