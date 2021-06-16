// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NAK_BSPLINE_BASE_DERIV2_HPP
#define NAK_BSPLINE_BASE_DERIV2_HPP

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
class NakBsplineBasisDeriv2 : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakBsplineBasisDeriv2() : bsplineBasis(BsplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineBasisDeriv2(size_t degree) : bsplineBasis(BsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakBsplineBasisDeriv2() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    double innerDeriv = hInvDbl * hInvDbl;
    double t = x * hInvDbl - static_cast<double>(i);

    switch (getDegree()) {
      case 1:
        return 0.0;

      case 3:
        if (l == 0) {
          return 0.0;
        } else if (l == 1) {
          if (i == 0) {
            // l = 0, i = 1
            return innerDeriv;
          } else if (i == 1) {
            // l = 1, i = 1
            return -2.0 * innerDeriv;
          } else {
            // l = 1, i = 2
            return innerDeriv;
          }
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 4, 3 < i < 2^l - 3
          return bsplineBasis.evalDxDx(l, i, x);
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
              double result = -2.5000000000000000e-01;
              result = 5.0000000000000000e-01 + result * t;
              return innerDeriv * result;
            }
          } else if ((l == 2) && (i == 1)) {
            // l = 2, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 5.9999999999999998e-01;
              result = -9.0000000000000002e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -1.4999999999999999e-01;
              result = 2.9999999999999999e-01 + result * t;
              return innerDeriv * result;
            }
          } else if (l == 2) {
            // l = 2, i = 2
            if ((t < -2.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 2.0;
              double result = -5.0000000000000000e-01;
              result = 4.0000000000000002e-01 + result * t;
              return innerDeriv * result;
            } else {
              double result = 5.0000000000000000e-01;
              result = -5.9999999999999998e-01 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 7.5000000000000000e-01;
              result = -1.0 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -5.0000000000000000e-01;
              result = 5.0000000000000000e-01 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 2) {
            // l >= 3, i = 2
            if ((t < -2.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 2.0;
              double result = -7.5000000000000000e-01;
              result = 5.0000000000000000e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 1.7500000000000000e+00;
              result = -1.0 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -7.5000000000000000e-01;
              result = 7.5000000000000000e-01 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 3, i = 3
            if ((t < -3.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 3.0;
              double result = 2.5000000000000000e-01;
              result *= t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -2.2500000000000000e+00;
              result = 5.0000000000000000e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 2.7500000000000000e+00;
              result = -1.7500000000000000e+00 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -1.0;
              result = 1.0 + result * t;
              return innerDeriv * result;
            }
          }
        }

      case 5:
        if (l == 0) {
          return 0.0;
        } else if (l == 1) {
          if (i == 0) {
            // l = 1, i = 0
            return innerDeriv;
          } else if (i == 1) {
            // l = 1, i = 1
            return -2.0 * innerDeriv;
          } else {
            // l = 1, i = 2
            return innerDeriv;
          }
        } else if ((i > 5) && (i < hInv - 5)) {
          // l >= 4, 5 < i < 2^l - 5
          return bsplineBasis.evalDxDx(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if ((l == 2) && (i == 0)) {
            // l = 2, i = 0
            double result = 1.0 / 2.0;
            result = -5.0 / 2.0 + result * t;
            result = 35.0 / 12.0 + result * t;
            return innerDeriv * result;
          } else if ((l == 2) && (i == 1)) {
            // l = 2, i = 1
            double result = -2.0;
            result = 5.0 + result * t;
            result = -5.0 / 3.0 + result * t;
            return innerDeriv * result;
          } else if ((l == 2) && (i == 2)) {
            // l = 2, i = 2
            double result = 3.0;
            result *= t;
            result = -5.0 / 2.0 + result * t;
            return innerDeriv * result;
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 7.0767195767195770e-02;
              result = -2.6984126984126983e-01 + result * t;
              result = 1.5873015873015872e-02 + result * t;
              result = 1.9576719576719576e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -2.6256613756613756e-01;
              result = 3.6706349206349204e-01 + result * t;
              result = 3.0753968253968256e-01 + result * t;
              result = -2.7447089947089948e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 1.5410052910052910e-01;
              result = -4.2063492063492064e-01 + result * t;
              result = 2.5396825396825395e-01 + result * t;
              result = 1.3756613756613756e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -4.6296296296296294e-03;
              result = 4.1666666666666664e-02 + result * t;
              result = -1.2500000000000000e-01 + result * t;
              result = 1.2500000000000000e-01 + result * t;
              return innerDeriv * result;
            }
          } else if ((l == 3) && (i == 4)) {
            // l = 3, i = 4
            if ((t < -4.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 4.0;
              double result = -2.7777777777777776e-02;
              result = 5.5555555555555552e-02 + result * t;
              result = 5.5555555555555552e-02 + result * t;
              result = 1.8518518518518517e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 2.5000000000000000e-01;
              result = -1.9444444444444445e-01 + result * t;
              result = -3.6111111111111110e-01 + result * t;
              result = -6.4814814814814811e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -2.5000000000000000e-01;
              result = 5.5555555555555558e-01 + result * t;
              result *= t;
              result = -3.7037037037037035e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = 2.7777777777777776e-02;
              result = -1.9444444444444445e-01 + result * t;
              result = 3.6111111111111110e-01 + result * t;
              result = -6.4814814814814811e-02 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 0) {
            // l >= 3, i = 0
            if ((t < 0.0) || (t > 3.0)) {
              return 0.0;
            } else {
              double result = -7.9365079365079361e-03;
              result = 7.1428571428571425e-02 + result * t;
              result = -2.1428571428571427e-01 + result * t;
              result = 2.1428571428571427e-01 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 3.9682539682539680e-02;
              result = -2.8571428571428570e-01 + result * t;
              result = 5.7142857142857140e-01 + result * t;
              result = -1.9047619047619047e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -2.3809523809523808e-02;
              result = 7.1428571428571425e-02 + result * t;
              result = -7.1428571428571425e-02 + result * t;
              result = 2.3809523809523808e-02 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 2) {
            // l >= 3, i = 2
            if ((t < -2.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 2.0;
              double result = -7.9365079365079361e-02;
              result = 4.2857142857142855e-01 + result * t;
              result = -4.2857142857142855e-01 + result * t;
              result = -2.3809523809523808e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 1.4285714285714285e-01;
              result = -2.8571428571428570e-01 + result * t;
              result *= t;
              result = 1.9047619047619047e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -4.7619047619047616e-02;
              result = 1.4285714285714285e-01 + result * t;
              result = -1.4285714285714285e-01 + result * t;
              result = 4.7619047619047616e-02 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 7.9365079365079361e-02;
              result = -2.8571428571428570e-01 + result * t;
              result *= t;
              result = 1.9047619047619047e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -3.6507936507936506e-01;
              result = 4.2857142857142855e-01 + result * t;
              result = 4.2857142857142855e-01 + result * t;
              result = -2.3809523809523808e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 3.0158730158730157e-01;
              result = -6.6666666666666663e-01 + result * t;
              result = 1.9047619047619047e-01 + result * t;
              result = 2.5396825396825395e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -7.9365079365079361e-02;
              result = 2.3809523809523808e-01 + result * t;
              result = -2.3809523809523808e-01 + result * t;
              result = 7.9365079365079361e-02 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 4) {
            // l >= 4, i = 4
            if ((t < -4.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 4.0;
              double result = -3.9682539682539680e-02;
              result = 7.1428571428571425e-02 + result * t;
              result = 7.1428571428571425e-02 + result * t;
              result = 2.3809523809523808e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 5.1587301587301593e-01;
              result = -2.8571428571428570e-01 + result * t;
              result = -5.7142857142857140e-01 + result * t;
              result = -1.9047619047619047e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -8.1746031746031744e-01;
              result = 1.2619047619047619e+00 + result * t;
              result = 4.0476190476190477e-01 + result * t;
              result = -5.3174603174603174e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 5.1587301587301593e-01;
              result = -1.1904761904761905e+00 + result * t;
              result = 4.7619047619047616e-01 + result * t;
              result = 3.1746031746031744e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -1.1904761904761904e-01;
              result = 3.5714285714285715e-01 + result * t;
              result = -3.5714285714285715e-01 + result * t;
              result = 1.1904761904761904e-01 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < -2.0) {
              t += 5.0;
              double result = 7.9365079365079361e-03;
              result *= t;
              result *= t;
              result *= t;
              return innerDeriv * result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = -4.3650793650793651e-01;
              result = 7.1428571428571425e-02 + result * t;
              result = 2.1428571428571427e-01 + result * t;
              result = 2.1428571428571427e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 1.2301587301587302e+00;
              result = -1.2380952380952381e+00 + result * t;
              result = -9.5238095238095233e-01 + result * t;
              result = 6.3492063492063489e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -1.4365079365079365e+00;
              result = 2.4523809523809526e+00 + result * t;
              result = 2.6190476190476192e-01 + result * t;
              result = -8.9682539682539686e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 7.8571428571428570e-01;
              result = -1.8571428571428572e+00 + result * t;
              result = 8.5714285714285710e-01 + result * t;
              result = 3.8095238095238093e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -1.6666666666666666e-01;
              result = 5.0000000000000000e-01 + result * t;
              result = -5.0000000000000000e-01 + result * t;
              result = 1.6666666666666666e-01 + result * t;
              return innerDeriv * result;
            }
          }
        }

      case 7:
        if (l == 0) {
          return 0.0;
        } else if (l == 1) {
          // l = 1, i = 1
          return -2.0 * innerDeriv;
        } else if ((i > 7) && (i < hInv - 7)) {
          // l >= 5, 5 < i < 2^l - 5
          return bsplineBasis.evalDxDx(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -2.0;
            result = 5.0 + result * t;
            result = -5.0 / 3.0 + result * t;
            return innerDeriv * result;
          } else if ((l == 3) && (i == 1)) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 2.3455710955710957e-04;
              result = -4.5211732711732712e-03 + result * t;
              result = 3.3449883449883450e-02 + result * t;
              result = -1.1204351204351204e-01 + result * t;
              result = 1.3706293706293707e-01 + result * t;
              result = 2.9588189588189588e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -8.4984459984459984e-06;
              result = 1.6996891996891997e-04 + result * t;
              result = -1.3597513597513598e-03 + result * t;
              result = 5.4390054390054390e-03 + result * t;
              result = -1.0878010878010878e-02 + result * t;
              result = 8.7024087024087024e-03 + result * t;
              return innerDeriv * result;
            }
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 5.1330613830613829e-04;
              result = -7.7165889665889668e-03 + result * t;
              result = 3.4537684537684535e-02 + result * t;
              result = -1.8492618492618493e-02 + result * t;
              result = -1.2346542346542347e-01 + result * t;
              result = 2.7521367521367520e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -2.1586052836052835e-04;
              result = 2.5495337995337995e-03 + result * t;
              result = -6.7987567987567990e-03 + result * t;
              result = -1.6317016317016316e-02 + result * t;
              result = 6.7987567987567984e-02 + result * t;
              result = -1.6317016317016317e-03 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 1) {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 4.8611111111111110e-04;
              result = -8.3333333333333332e-03 + result * t;
              result = 5.2777777777777778e-02 + result * t;
              result = -1.4166666666666666e-01 + result * t;
              result = 1.0972222222222222e-01 + result * t;
              result = 8.5833333333333331e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -2.7777777777777778e-04;
              result = 1.3888888888888889e-03 + result * t;
              result = -2.7777777777777779e-03 + result * t;
              result = 2.7777777777777779e-03 + result * t;
              result = -1.3888888888888889e-03 + result * t;
              result = 2.7777777777777778e-04 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 3.0;
              double result = 2.4305555555555556e-03;
              result = -2.7777777777777776e-02 + result * t;
              result = 8.3333333333333329e-02 + result * t;
              result = 2.7777777777777776e-02 + result * t;
              result = -2.5000000000000000e-01 + result * t;
              result = -4.7222222222222221e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -9.0277777777777769e-03;
              result = 2.0833333333333332e-02 + result * t;
              result = 2.7777777777777776e-02 + result * t;
              result = -8.3333333333333329e-02 + result * t;
              result = -2.7777777777777776e-02 + result * t;
              result = 1.0833333333333334e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 6.2500000000000003e-03;
              result = -2.4305555555555556e-02 + result * t;
              result = 2.0833333333333332e-02 + result * t;
              result = 3.4722222222222224e-02 + result * t;
              result = -7.2916666666666671e-02 + result * t;
              result = 3.6805555555555557e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -1.3888888888888889e-03;
              result = 6.9444444444444441e-03 + result * t;
              result = -1.3888888888888888e-02 + result * t;
              result = 1.3888888888888888e-02 + result * t;
              result = -6.9444444444444441e-03 + result * t;
              result = 1.3888888888888889e-03 + result * t;
              return innerDeriv * result;
            }
          } else if (i == 5) {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 5.0;
              double result = 1.4583333333333334e-03;
              result = -8.3333333333333332e-03 + result * t;
              result = -2.7777777777777779e-03 + result * t;
              result = 2.5000000000000001e-02 + result * t;
              result = 4.0277777777777780e-02 + result * t;
              result = 1.9166666666666665e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -3.0624999999999999e-02;
              result = 2.0833333333333332e-02 + result * t;
              result = 9.7222222222222224e-02 + result * t;
              result = 1.2500000000000000e-01 + result * t;
              result = -1.5972222222222221e-01 + result * t;
              result = -2.3749999999999999e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 6.1041666666666668e-02;
              result = -1.3229166666666667e-01 + result * t;
              result = -1.2569444444444444e-01 + result * t;
              result = 2.3541666666666666e-01 + result * t;
              result = 3.1215277777777778e-01 + result * t;
              result = -1.8479166666666666e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -5.3541666666666668e-02;
              result = 1.7291666666666666e-01 + result * t;
              result = -4.4444444444444446e-02 + result * t;
              result = -3.2500000000000001e-01 + result * t;
              result = 1.8194444444444444e-01 + result * t;
              result = 1.6583333333333333e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 2.2847222222222224e-02;
              result = -9.4791666666666663e-02 + result * t;
              result = 1.1180555555555556e-01 + result * t;
              result = 4.3749999999999997e-02 + result * t;
              result = -1.7743055555555556e-01 + result * t;
              result = 9.7708333333333328e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -3.8888888888888888e-03;
              result = 1.9444444444444445e-02 + result * t;
              result = -3.8888888888888890e-02 + result * t;
              result = 3.8888888888888890e-02 + result * t;
              result = -1.9444444444444445e-02 + result * t;
              result = 3.8888888888888888e-03 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 7
            if ((t < -7.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < -3.0) {
              t += 7.0;
              double result = 6.9444444444444444e-05;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              return innerDeriv * result;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -2.2847222222222224e-02;
              result = 1.3888888888888889e-03 + result * t;
              result = 1.1111111111111112e-02 + result * t;
              result = 4.4444444444444446e-02 + result * t;
              result = 8.8888888888888892e-02 + result * t;
              result = 7.1111111111111111e-02 + result * t;
              return innerDeriv * result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 1.0548611111111111e-01;
              result = -1.1284722222222222e-01 + result * t;
              result = -2.1180555555555555e-01 + result * t;
              result = -1.4236111111111110e-01 + result * t;
              result = 1.0243055555555555e-01 + result * t;
              result = 1.9409722222222223e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -2.1534722222222222e-01;
              result = 4.1458333333333336e-01 + result * t;
              result = 3.9166666666666666e-01 + result * t;
              result = -4.0000000000000002e-01 + result * t;
              result = -7.4166666666666670e-01 + result * t;
              result = -6.5000000000000002e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = 2.4298611111111112e-01;
              result = -6.6215277777777781e-01 + result * t;
              result = -1.0347222222222222e-01 + result * t;
              result = 1.1090277777777777e+00 + result * t;
              result = 2.1493055555555557e-01 + result * t;
              result = -6.1576388888888889e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -1.5805555555555556e-01;
              result = 5.5277777777777781e-01 + result * t;
              result = -3.2222222222222224e-01 + result * t;
              result = -7.4444444444444446e-01 + result * t;
              result = 6.8888888888888888e-01 + result * t;
              result = 1.8555555555555556e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 5.5833333333333332e-02;
              result = -2.3749999999999999e-01 + result * t;
              result = 3.0833333333333335e-01 + result * t;
              result = 2.5000000000000001e-02 + result * t;
              result = -3.4583333333333333e-01 + result * t;
              result = 2.0250000000000001e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -8.3333333333333332e-03;
              result = 4.1666666666666664e-02 + result * t;
              result = -8.3333333333333329e-02 + result * t;
              result = 8.3333333333333329e-02 + result * t;
              result = -4.1666666666666664e-02 + result * t;
              result = 8.3333333333333332e-03 + result * t;
              return innerDeriv * result;
            }
          }
        }

      default:
        return 0.0;
    }
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "NakBsplineBasisDeriv2: evalDx not implemented" << std::endl;
    return -1;
  }

  inline double getIntegral(LT level, IT index) override {
    std::cerr << "NakBsplineBasisDeriv2: getIntegral not implemented" << std::endl;
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
typedef NakBsplineBasisDeriv2<unsigned int, unsigned int> SNakBsplineBaseDeriv2;

}  // namespace base
}  // namespace sgpp

#endif /* NAK_BSPLINE_BASE_DERIV2_HPP */
