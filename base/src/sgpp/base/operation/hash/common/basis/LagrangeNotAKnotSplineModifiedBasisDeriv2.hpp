// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LAGRANGE_NOTAKNOT_SPLINE_MODIFIED_BASE_DERIV2_HPP
#define LAGRANGE_NOTAKNOT_SPLINE_MODIFIED_BASE_DERIV2_HPP

#include <sgpp/base/operation/hash/common/basis/LagrangeNotAKnotSplineBasisDeriv2.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <cstddef>
#include <algorithm>

namespace sgpp {
namespace base {

/**
 * Lagrange not-a-knot spline basis.
 */
template <class LT, class IT>
class LagrangeNotAKnotSplineModifiedBasisDeriv2: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  LagrangeNotAKnotSplineModifiedBasisDeriv2() :
    lagrangeNotAKnotSplineBasisDeriv2(LagrangeNotAKnotSplineBasisDeriv2<LT, IT>()) {
  }

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit LagrangeNotAKnotSplineModifiedBasisDeriv2(size_t degree) :
    lagrangeNotAKnotSplineBasisDeriv2(LagrangeNotAKnotSplineBasisDeriv2<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~LagrangeNotAKnotSplineModifiedBasisDeriv2() override {
  }

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
        return lagrangeNotAKnotSplineBasisDeriv2.eval(l, i, x);

      case 3:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return lagrangeNotAKnotSplineBasisDeriv2.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 2.5961538461538464e-01;
              result *= t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -3.4615384615384615e-01;
              result = 5.1923076923076927e-01 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 3.2142857142857145e-01;
              result *= t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -7.5000000000000000e-01;
              result = 6.4285714285714290e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = 1.0714285714285714e-01;
              result = -1.0714285714285714e-01 + result * t;
              return innerDeriv * result;
            }
          }
        }

      case 5:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return lagrangeNotAKnotSplineBasisDeriv2.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -18.0/35.0;
            result = 19.0/35.0 + result * t;
            result = 37.0/35.0 + result * t;
            return innerDeriv * result;
          } else if (l == 3) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 4.3883239512345067e-02;
              result = -2.8687446732898053e-01 + result * t;
              result = 4.6545741193536394e-01 + result * t;
              result *= t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -4.2262363968840908e-02;
              result = 1.0807468828212510e-01 + result * t;
              result = -7.0941925205202289e-02 + result * t;
              result = -6.5050332141593551e-04 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 6.2632990450480746e-03;
              result = -1.8712403624397633e-02 + result * t;
              result = 1.8420359452525181e-02 + result * t;
              result = -5.7801042133340313e-03 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -9.1168836172461064e-06;
              result = 7.7493510746591896e-05 + result * t;
              result = -2.1455066112585836e-04 + result * t;
              result = 1.9115065984159335e-04 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 4.4192004314703701e-02;
              result = -2.8841097552753997e-01 + result * t;
              result = 4.6706667565563270e-01 + result * t;
              result *= t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -4.3180745634733363e-02;
              result = 1.0931706330479336e-01 + result * t;
              result = -7.0215061012607022e-02 + result * t;
              result = -1.3146362839614374e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 6.9664486842400952e-03;
              result = -2.0225173599406728e-02 + result * t;
              result = 1.8876828692779613e-02 + result * t;
              result = -5.3933796265084608e-03 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -2.2472415110451921e-04;
              result = 6.7417245331355760e-04 + result * t;
              result = -6.7417245331355760e-04 + result * t;
              result = 2.2472415110451921e-04 + result * t;
              return innerDeriv * result;
            }
          }
        }

      case 7:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return lagrangeNotAKnotSplineBasisDeriv2.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -18.0/35.0;
            result = 19.0/35.0 + result * t;
            result = 37.0/35.0 + result * t;
            return innerDeriv * result;
          } else if (l == 3) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 4.9273252509971508e-04;
              result = -9.1337602465057875e-03 + result * t;
              result = 6.3797769812555707e-02 + result * t;
              result = -1.9795489764402116e-01 + result * t;
              result = 2.2836666498747241e-01 + result * t;
              result *= t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -5.0790597792481687e-05;
              result = 7.2089025548851269e-04 + result * t;
              result = -3.5051901155824920e-03 + result * t;
              result = 6.1261725059093414e-03 + result * t;
              result = -5.2455613986926691e-04 + result * t;
              result = -4.4389517542570622e-03 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 7.6455663668362130e-04;
              result = -1.2901695723965144e-02 + result * t;
              result = 8.1144297709566521e-02 + result * t;
              result = -2.2480107290108686e-01 + result * t;
              result = 2.3054729004069041e-01 + result * t;
              result *= t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -5.7553345725583035e-04;
              result = 2.3894370097072805e-03 + result * t;
              result = -2.9537720044963910e-03 + result * t;
              result = -3.1604240942486824e-04 + result * t;
              result = 2.8633865111468005e-03 + result * t;
              result = -1.3210622134195524e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 1.0408813751938767e-04;
              result = -4.8823027657187107e-04 + result * t;
              result = 8.4864146177442764e-04 + result * t;
              result = -5.9607093722866206e-04 + result * t;
              result = 5.0066431357860826e-05 + result * t;
              result = 8.6413436257438717e-05 + result * t;
              return innerDeriv * result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -6.4932098415611783e-06;
              result = 3.2210411025067265e-05 + result * t;
              result = -6.3398269319180013e-05 + result * t;
              result = 6.1353163857270980e-05 + result * t;
              result = -2.8631476466726456e-05 + result * t;
              result = 4.9082531085816785e-06 + result * t;
              return innerDeriv * result;
            } else {
              t -= 6.0;
              double result = 5.1127636547725818e-08;
              result = -2.5563818273862909e-07 + result * t;
              result = 5.1127636547725818e-07 + result * t;
              result = -5.1127636547725818e-07 + result * t;
              result = 2.5563818273862909e-07 + result * t;
              result = -5.1127636547725818e-08 + result * t;
              return innerDeriv * result;
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
    return lagrangeNotAKnotSplineBasisDeriv2.getDegree();
  }

 protected:
  /// Unmodified basis
  LagrangeNotAKnotSplineBasisDeriv2<LT, IT> lagrangeNotAKnotSplineBasisDeriv2;
};

// default type-def (unsigned int for level and index)
typedef LagrangeNotAKnotSplineModifiedBasisDeriv2<unsigned int, unsigned int>
SLagrangeNotAKnotSplineModifiedBaseDeriv2;

}  // namespace base
}  // namespace sgpp

#endif /* LAGRANGE_NOTAKNOT_SPLINE_MODIFIED_BASE_DERIV2_HPP */
