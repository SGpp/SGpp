// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NAK_BSPLINE_MODIFIED_BASE_DERIV1_HPP
#define NAK_BSPLINE_MODIFIED_BASE_DERIV1_HPP

#include <sgpp/base/operation/hash/common/basis/NakBsplineBasisDeriv1.hpp>

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
            double result = -6.0 / 35.0;
            result = 19.0 / 70.0 + result * t;
            result = 37.0 / 35.0 + result * t;
            result = -331.0 / 210.0 + result * t;
            return innerDeriv * result;
          } else {
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
          }
        }

      case 7:
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
            double result = -6.0 / 35.0;
            result = 19.0 / 70.0 + result * t;
            result = 37.0 / 35.0 + result * t;
            result = -331.0 / 210.0 + result * t;
            return innerDeriv * result;
          } else if (l == 3) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 4.3908637658637657e-05;
              result = -1.0198135198135198e-03 + result * t;
              result = 9.5182595182595180e-03 + result * t;
              result = -4.3512043512043512e-02 + result * t;
              result = 8.7024087024087024e-02 + result * t;
              result *= t;
              result = -1.8565138565138564e-01 + result * t;
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
          } else {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 9.4988787615740740e-05;
              result = -2.0019531249999998e-03 + result * t;
              result = 1.6547309027777776e-02 + result * t;
              result = -6.5104166666666671e-02 + result * t;
              result = 1.0850694444444445e-01 + result * t;
              result *= t;
              result = -1.4467592592592593e-01 + result * t;
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
          }
        }

      default:
        return 0.0;
    }
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "NakBsplineModifiedBasisDeriv1: evalDx not implemented" << std::endl;
    return -1;
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

#endif /* NAK_BSPLINE_MODIFIED_BASE_DERIV1_HPP */
