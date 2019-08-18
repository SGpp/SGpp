// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NAK_BSPLINE_MODIFIED_BASE_DERIV2_HPP
#define NAK_BSPLINE_MODIFIED_BASE_DERIV2_HPP

#include <sgpp/base/operation/hash/common/basis/NakBsplineBasisDeriv2.hpp>

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
class NakBsplineModifiedBasisDeriv2: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakBsplineModifiedBasisDeriv2() :
    nakBsplineBasisDeriv2(NakBsplineBasisDeriv2<LT, IT>()) {
  }

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineModifiedBasisDeriv2(size_t degree) :
        nakBsplineBasisDeriv2(NakBsplineBasisDeriv2<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakBsplineModifiedBasisDeriv2() override {
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
        return 0.0;

      case 3:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return nakBsplineBasisDeriv2.eval(l, i, x);
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
              double result = 1.4999999999999999e-01;
              result *= t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -1.4999999999999999e-01;
              result = 2.9999999999999999e-01 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 2.5000000000000000e-01;
              result *= t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -5.0000000000000000e-01;
              result = 5.0000000000000000e-01 + result * t;
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
          return nakBsplineBasisDeriv2.eval(l, i, x);
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
          } else {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 3.2627865961199293e-02;
              result = -2.2222222222222221e-01 + result * t;
              result = 3.8095238095238093e-01 + result * t;
              result *= t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -2.3809523809523808e-02;
              result = 7.1428571428571425e-02 + result * t;
              result = -7.1428571428571425e-02 + result * t;
              result = 2.3809523809523808e-02 + result * t;
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
          return nakBsplineBasisDeriv2.eval(l, i, x);
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
              double result = 2.6345182595182597e-04;
              result = -5.0990675990675990e-03 + result * t;
              result = 3.8073038073038072e-02 + result * t;
              result = -1.3053613053613053e-01 + result * t;
              result = 1.7404817404817405e-01 + result * t;
              result *= t;
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
          } else {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 5.6993272569444441e-04;
              result = -1.0009765625000000e-02 + result * t;
              result = 6.6189236111111105e-02 + result * t;
              result = -1.9531250000000000e-01 + result * t;
              result = 2.1701388888888890e-01 + result * t;
              result *= t;
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
    return nakBsplineBasisDeriv2.getDegree();
  }

 protected:
  /// B-spline basis for B-spline evaluation
  NakBsplineBasisDeriv2<LT, IT> nakBsplineBasisDeriv2;
};

// default type-def (unsigned int for level and index)
typedef NakBsplineModifiedBasisDeriv2<unsigned int, unsigned int>
SNakBsplineModifiedBaseDeriv2;

}  // namespace base
}  // namespace sgpp

#endif /* NAK_BSPLINE_MODIFIED_BASE_DERIV2_HPP */
