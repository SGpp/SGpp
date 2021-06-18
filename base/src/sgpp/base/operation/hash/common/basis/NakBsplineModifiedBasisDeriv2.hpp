// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NAK_BSPLINE_MODIFIED_BASE_DERIV2_HPP
#define NAK_BSPLINE_MODIFIED_BASE_DERIV2_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasisDeriv2.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Not-a-knot B-spline basis.
 */
template <class LT, class IT>
class NakBsplineModifiedBasisDeriv2 : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakBsplineModifiedBasisDeriv2() : nakBsplineBasisDeriv2(NakBsplineBasisDeriv2<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineModifiedBasisDeriv2(size_t degree)
      : nakBsplineBasisDeriv2(NakBsplineBasisDeriv2<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakBsplineModifiedBasisDeriv2() override {}

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
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 3, 3 < i < 2^l - 3
          return nakBsplineBasisDeriv2.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1 (Lagrange)
            double result = -6.0 / 23.0;
            result = 12.0 / 23.0 + result * t;
            result = 18.0 / 23.0 + result * t;
            return innerDeriv * result;
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i=3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 5.5114638447971778e-03;
              result = 8.2539682539682538e-02 + result * t;
              result = -3.3650793650793653e-01 + result * t;
              result *= t;
              return innerDeriv * result;
            } else if (t < 1.0) {
              double result = -1.4510582010582010e-01;
              result = 1.3214285714285715e-01 + result * t;
              result = 3.0753968253968256e-01 + result * t;
              result = -1.1785714285714285e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 1.1494708994708995e-01;
              result = -3.0317460317460315e-01 + result * t;
              result = 1.3650793650793649e-01 + result * t;
              result = 1.7671957671957672e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -4.6296296296296294e-03;
              result = 4.1666666666666664e-02 + result * t;
              result = -1.2500000000000000e-01 + result * t;
              result = 1.2500000000000000e-01 + result * t;
              return innerDeriv * result;
            }
          } else {
            if (i == 1) {
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
            } else {
              // l >= 4, i = 3
              if ((t < -3.0) || (t > 3.0)) {
                return 0.0;
              } else if (t < 0.0) {
                t += 3.0;
                double result = 1.0 / 63.0;
                result = 2.0 / 35.0 + result * t;
                result = -12.0 / 35.0 + result * t;
                result *= t;
                return innerDeriv * result;
              } else if (t < 1.0) {
                double result = -79.0 / 315.0;
                result = 1.0 / 5.0 + result * t;
                result = 3.0 / 7.0 + result * t;
                result = -3.0 / 35.0 + result * t;
                return innerDeriv * result;
              } else if (t < 2.0) {
                t -= 1.0;
                double result = 83.0 / 315.0;
                result = -58.0 / 105.0 + result * t;
                result = 8.0 / 105.0 + result * t;
                result = 92.0 / 315.0 + result * t;
                return innerDeriv * result;
              } else {
                t -= 2.0;
                double result = -5.0 / 63.0;
                result = 5.0 / 21.0 + result * t;
                result = -5.0 / 21.0 + result * t;
                result = 5.0 / 63.0 + result * t;
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
    std::cerr << "NakBsplineMOdifiedBasisDeriv2::evalDx not implemented\n";
    return 0;
  }

  inline double getIntegral(LT level, IT index) override { return -1.0; }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override { return nakBsplineBasisDeriv2.getDegree(); }

 protected:
  /// B-spline basis for B-spline evaluation
  NakBsplineBasisDeriv2<LT, IT> nakBsplineBasisDeriv2;
};

// default type-def (unsigned int for level and index)
typedef NakBsplineModifiedBasisDeriv2<unsigned int, unsigned int> SNakBsplineModifiedBaseDeriv2;

}  // namespace base
}  // namespace sgpp

#endif /* NAK_BSPLINE_MODIFIED_BASE_DERIV2_HPP */
