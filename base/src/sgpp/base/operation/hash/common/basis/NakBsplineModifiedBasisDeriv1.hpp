// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasisDeriv1.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Not-a-knot B-spline basis.
 */
template <class LT, class IT>
class NakBsplineModifiedBasisDeriv1 : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakBsplineModifiedBasisDeriv1() : nakBsplineBasisDeriv1(NakBsplineBasisDeriv1<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineModifiedBasisDeriv1(size_t degree)
      : nakBsplineBasisDeriv1(NakBsplineBasisDeriv1<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakBsplineModifiedBasisDeriv1() override {}

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
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return nakBsplineBasisDeriv1.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }

          // l >= 3, i = 1
          if (t > 1.0) {
            return 0.0;
          } else {
            return -innerDeriv;
          }
        }

      case 3:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return nakBsplineBasisDeriv1.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 7.4999999999999997e-02;
              result *= t;
              result = -5.9999999999999998e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -7.4999999999999997e-02;
              result = 2.9999999999999999e-01 + result * t;
              result = -2.9999999999999999e-01 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 1.2500000000000000e-01;
              result *= t;
              result = -7.5000000000000000e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -2.5000000000000000e-01;
              result = 5.0000000000000000e-01 + result * t;
              result = -2.5000000000000000e-01 + result * t;
              return innerDeriv * result;
            }
          }
        }

      case 5:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 3, 3 < i < 2^l - 3
          return nakBsplineBasisDeriv1.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1 (Lagrange)
            double result = -2.0 / 23.0;
            result = 6.0 / 23.0 + result * t;
            result = 18.0 / 23.0 + result * t;
            result = -67.0 / 46.0 + result * t;
            return innerDeriv * result;
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i=3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 1.3778659611992945e-03;
              result = 2.7513227513227514e-02 + result * t;
              result = -1.6825396825396827e-01 + result * t;
              result *= t;
              result = 3.4973544973544973e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -3.6276455026455025e-02;
              result = 4.4047619047619051e-02 + result * t;
              result = 1.5376984126984128e-01 + result * t;
              result = -1.1785714285714285e-01 + result * t;
              result = -3.1008597883597883e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 2.8736772486772488e-02;
              result = -1.0105820105820106e-01 + result * t;
              result = 6.8253968253968247e-02 + result * t;
              result = 1.7671957671957672e-01 + result * t;
              result = -2.6640211640211642e-01 + result * t;
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
          } else {
            if (i == 1) {
              // l >= 3, i = 1
              if ((t < -1.0) || (t > 3.0)) {
                return 0.0;
              } else if (t < 2.0) {
                t += 1.0;
                double result = 8.1569664902998232e-03;
                result = -7.4074074074074070e-02 + result * t;
                result = 1.9047619047619047e-01 + result * t;
                result *= t;
                result = -3.8095238095238093e-01 + result * t;
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
            } else {
              // l >= 4, i = 3
              if ((t < -3.0) || (t > 3.0)) {
                return 0.0;
              } else if (t < 0.0) {
                t += 3.0;
                double result = 3.9682539682539680e-03;
                result = 1.9047619047619049e-02 + result * t;
                result = -1.7142857142857143e-01 + result * t;
                result *= t;
                result = 3.4285714285714286e-01 + result * t;
                return innerDeriv * result;
              } else if (t < 1.0) {
                double result = -6.2698412698412698e-02;
                result = 6.6666666666666666e-02 + result * t;
                result = 2.1428571428571427e-01 + result * t;
                result = -8.5714285714285715e-02 + result * t;
                result = -3.6428571428571427e-01 + result * t;
                return innerDeriv * result;
              } else if (t < 2.0) {
                t -= 1.0;
                double result = 6.5873015873015875e-02;
                result = -1.8412698412698414e-01 + result * t;
                result = 3.8095238095238099e-02 + result * t;
                result = 2.9206349206349208e-01 + result * t;
                result = -2.3174603174603176e-01 + result * t;
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
   * @return      value of derivative of wavelet basis function
   */
  inline double evalDx(LT l, IT i, double x) override {
    std::cerr << "NakBsplineMOdifiedBasisDeriv1::evalDx not implemented\n";
    return 0;
  }

  inline double getIntegral(LT level, IT index) override {
    std::cerr << "NakBsplineModifiedBasisDeriv1: getIntegral not implemented" << std::endl;
    return -1.0;
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override { return nakBsplineBasisDeriv1.getDegree(); }

 protected:
  /// B-spline basis for B-spline evaluation
  NakBsplineBasisDeriv1<LT, IT> nakBsplineBasisDeriv1;
};

// default type-def (unsigned int for level and index)
typedef NakBsplineModifiedBasisDeriv1<unsigned int, unsigned int> SNakBsplineModifiedBaseDeriv1;

}  // namespace base
}  // namespace sgpp
