// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NAK_BSPLINE_BASE_DERIV1_HPP
#define NAK_BSPLINE_BASE_DERIV1_HPP

#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * Not-a-knot B-spline basis.
 */
template <class LT, class IT>
class NakBsplineBasisDeriv1 : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakBsplineBasisDeriv1() : bsplineBasis(BsplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineBasisDeriv1(size_t degree) : bsplineBasis(BsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakBsplineBasisDeriv1() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double eval(LT l, IT i, double x) override {
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

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "NakBsplineBasisDeriv1: evalDx not implemented" << std::endl;
    return -1;
  }

  inline double getIntegral(LT level, IT index) override {
    std::cerr << "NakBsplineBasisDeriv1: getIntegral not implemented" << std::endl;
    return -1.0;
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override { return bsplineBasis.getDegree(); }

 protected:
  /// B-spline basis for B-spline evaluation
  BsplineBasis<LT, IT> bsplineBasis;
};

// default type-def (unsigned int for level and index)
typedef NakBsplineBasisDeriv1<unsigned int, unsigned int> SNakBsplineBaseDeriv1;

}  // namespace base
}  // namespace sgpp

#endif /* NAK_BSPLINE_BASE_DERIV1_HPP */
