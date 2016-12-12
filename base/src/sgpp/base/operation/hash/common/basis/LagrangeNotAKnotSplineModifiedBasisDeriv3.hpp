// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LAGRANGE_NOTAKNOT_SPLINE_MODIFIED_BASE_DERIV3_HPP
#define LAGRANGE_NOTAKNOT_SPLINE_MODIFIED_BASE_DERIV3_HPP

#include <sgpp/base/operation/hash/common/basis/LagrangeNotAKnotSplineBasisDeriv3.hpp>

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
class LagrangeNotAKnotSplineModifiedBasisDeriv3: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  LagrangeNotAKnotSplineModifiedBasisDeriv3() :
    lagrangeNotAKnotSplineBasisDeriv3(LagrangeNotAKnotSplineBasisDeriv3<LT, IT>()) {
  }

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit LagrangeNotAKnotSplineModifiedBasisDeriv3(size_t degree) :
    lagrangeNotAKnotSplineBasisDeriv3(LagrangeNotAKnotSplineBasisDeriv3<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~LagrangeNotAKnotSplineModifiedBasisDeriv3() override {
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
    double innerDeriv = hInvDbl * hInvDbl * hInvDbl;
    double t = x * hInvDbl - static_cast<double>(i);

    switch (getDegree()) {
      case 1:
        return lagrangeNotAKnotSplineBasisDeriv3.eval(l, i, x);

      case 3:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return lagrangeNotAKnotSplineBasisDeriv3.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= 1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              return innerDeriv * 2.5961538461538464e-01;
            } else {
              return innerDeriv * -3.4615384615384615e-01;
            }
          } else {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              return innerDeriv * 3.2142857142857145e-01;
            } else if (t < 2.0) {
              return innerDeriv * -7.5000000000000000e-01;
            } else {
              return innerDeriv * 1.0714285714285714e-01;
            }
          }
        }

      case 5:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return lagrangeNotAKnotSplineBasisDeriv3.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= 1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -36.0/35.0;
            result = 19.0/35.0 + result * t;
            return innerDeriv * result;
          } else if (l == 3) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 1.3164971853703519e-01;
              result = -5.7374893465796106e-01 + result * t;
              result = 4.6545741193536394e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -1.2678709190652274e-01;
              result = 2.1614937656425021e-01 + result * t;
              result = -7.0941925205202289e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 1.8789897135144226e-02;
              result = -3.7424807248795267e-02 + result * t;
              result = 1.8420359452525181e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -2.7350650851738319e-05;
              result = 1.5498702149318379e-04 + result * t;
              result = -2.1455066112585836e-04 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 1.3257601294411112e-01;
              result = -5.7682195105507994e-01 + result * t;
              result = 4.6706667565563270e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -1.2954223690420008e-01;
              result = 2.1863412660958673e-01 + result * t;
              result = -7.0215061012607022e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 2.0899346052720286e-02;
              result = -4.0450347198813456e-02 + result * t;
              result = 1.8876828692779613e-02 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -6.7417245331355760e-04;
              result = 1.3483449066271152e-03 + result * t;
              result = -6.7417245331355760e-04 + result * t;
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
          return lagrangeNotAKnotSplineBasisDeriv3.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= 1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -36.0/35.0;
            result = 19.0/35.0 + result * t;
            return innerDeriv * result;
          } else if (l == 3) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 2.4636626254985754e-03;
              result = -3.6535040986023150e-02 + result * t;
              result = 1.9139330943766714e-01 + result * t;
              result = -3.9590979528804232e-01 + result * t;
              result = 2.2836666498747241e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -2.5395298896240844e-04;
              result = 2.8835610219540508e-03 + result * t;
              result = -1.0515570346747476e-02 + result * t;
              result = 1.2252345011818683e-02 + result * t;
              result = -5.2455613986926691e-04 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 3.8227831834181064e-03;
              result = -5.1606782895860577e-02 + result * t;
              result = 2.4343289312869956e-01 + result * t;
              result = -4.4960214580217372e-01 + result * t;
              result = 2.3054729004069041e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -2.8776672862791514e-03;
              result = 9.5577480388291218e-03 + result * t;
              result = -8.8613160134891726e-03 + result * t;
              result = -6.3208481884973649e-04 + result * t;
              result = 2.8633865111468005e-03 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 5.2044068759693831e-04;
              result = -1.9529211062874843e-03 + result * t;
              result = 2.5459243853232826e-03 + result * t;
              result = -1.1921418744573241e-03 + result * t;
              result = 5.0066431357860826e-05 + result * t;
              return innerDeriv * result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -3.2466049207805891e-05;
              result = 1.2884164410026906e-04 + result * t;
              result = -1.9019480795754004e-04 + result * t;
              result = 1.2270632771454196e-04 + result * t;
              result = -2.8631476466726456e-05 + result * t;
              return innerDeriv * result;
            } else {
              t -= 6.0;
              double result = 2.5563818273862909e-07;
              result = -1.0225527309545164e-06 + result * t;
              result = 1.5338290964317745e-06 + result * t;
              result = -1.0225527309545164e-06 + result * t;
              result = 2.5563818273862909e-07 + result * t;
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
    return lagrangeNotAKnotSplineBasisDeriv3.getDegree();
  }

 protected:
  /// Unmodified basis
  LagrangeNotAKnotSplineBasisDeriv3<LT, IT> lagrangeNotAKnotSplineBasisDeriv3;
};

// default type-def (unsigned int for level and index)
typedef LagrangeNotAKnotSplineModifiedBasisDeriv3<unsigned int, unsigned int>
SLagrangeNotAKnotSplineModifiedBaseDeriv3;

}  // namespace base
}  // namespace sgpp

#endif /* LAGRANGE_NOTAKNOT_SPLINE_MODIFIED_BASE_DERIV3_HPP */
