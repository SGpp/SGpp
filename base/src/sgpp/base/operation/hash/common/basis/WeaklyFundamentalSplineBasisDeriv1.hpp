// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WEAKLY_FUNDAMENTAL_SPLINE_BASE_DERIV1_HPP
#define WEAKLY_FUNDAMENTAL_SPLINE_BASE_DERIV1_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * weakly fundamental spline basis (1st derivative).
 */
template <class LT, class IT>
class WeaklyFundamentalSplineBasisDeriv1 : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  WeaklyFundamentalSplineBasisDeriv1() : degree(0) {}

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit WeaklyFundamentalSplineBasisDeriv1(size_t degree) : degree(degree) {
    if (degree < 1) {
      this->degree = 1;
    } else if (degree % 2 == 0) {
      this->degree = degree - 1;
    }

    if (this->degree > 7) {
      throw std::runtime_error("Unsupported weakly fundamental spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~WeaklyFundamentalSplineBasisDeriv1() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    const double innerDeriv = hInvDbl;
    double t = x * hInvDbl - static_cast<double>(i);

    switch (degree) {
      case 1:
        if ((t < -1.0) || (t > 1.0)) {
          return 0.0;
        } else if (t < 0.0) {
          return innerDeriv;
        } else {
          return -innerDeriv;
        }

      case 3:
        if ((t < -3.0) || (t > 3.0)) {
          return 0.0;
        } else if (t < -2.0) {
          t += 3.0;
          double result = -1.2500000000000000e-01;
          result *= t;
          result *= t;
          return innerDeriv * result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = 8.7500000000000000e-01;
          result = -2.5000000000000000e-01 + result * t;
          result = -1.2500000000000000e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = -2.0;
          result = 1.5000000000000000e+00 + result * t;
          result = 5.0000000000000000e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 1.0) {
          double result = 2.0;
          result = -2.5000000000000000e+00 + result * t;
          result *= t;
          return innerDeriv * result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = -8.7500000000000000e-01;
          result = 1.5000000000000000e+00 + result * t;
          result = -5.0000000000000000e-01 + result * t;
          return innerDeriv * result;
        } else {
          t -= 2.0;
          double result = 1.2500000000000000e-01;
          result = -2.5000000000000000e-01 + result * t;
          result = 1.2500000000000000e-01 + result * t;
          return innerDeriv * result;
        }

      case 5:
        if ((t < -5.0) || (t > 5.0)) {
          return 0.0;
        } else if (t < -4.0) {
          t += 5.0;
          double result = 6.3131313131313137e-04;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          return innerDeriv * result;
        } else if (t < -3.0) {
          t += 4.0;
          double result = -1.9570707070707072e-02;
          result = 2.5252525252525255e-03 + result * t;
          result = 3.7878787878787880e-03 + result * t;
          result = 2.5252525252525255e-03 + result * t;
          result = 6.3131313131313137e-04 + result * t;
          return innerDeriv * result;
        } else if (t < -2.0) {
          t += 3.0;
          double result = 1.3005050505050506e-01;
          result = -7.5757575757575760e-02 + result * t;
          result = -1.0606060606060606e-01 + result * t;
          result = -6.0606060606060608e-02 + result * t;
          result = -1.0101010101010102e-02 + result * t;
          return innerDeriv * result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = -3.9520202020202022e-01;
          result = 4.4444444444444442e-01 + result * t;
          result = 4.4696969696969696e-01 + result * t;
          result = 2.0202020202020204e-02 + result * t;
          result = -1.2247474747474747e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = 6.6666666666666663e-01;
          result = -1.1363636363636365e+00 + result * t;
          result = -5.9090909090909094e-01 + result * t;
          result = 6.6666666666666663e-01 + result * t;
          result = 3.9393939393939392e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 1.0) {
          double result = -6.6666666666666663e-01;
          result = 1.5303030303030303e+00 + result * t;
          result *= t;
          result = -1.2575757575757576e+00 + result * t;
          result *= t;
          return innerDeriv * result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = 3.9520202020202022e-01;
          result = -1.1363636363636365e+00 + result * t;
          result = 5.9090909090909094e-01 + result * t;
          result = 6.6666666666666663e-01 + result * t;
          result = -3.9393939393939392e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 3.0) {
          t -= 2.0;
          double result = -1.3005050505050506e-01;
          result = 4.4444444444444442e-01 + result * t;
          result = -4.4696969696969696e-01 + result * t;
          result = 2.0202020202020204e-02 + result * t;
          result = 1.2247474747474747e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 4.0) {
          t -= 3.0;
          double result = 1.9570707070707072e-02;
          result = -7.5757575757575760e-02 + result * t;
          result = 1.0606060606060606e-01 + result * t;
          result = -6.0606060606060608e-02 + result * t;
          result = 1.0101010101010102e-02 + result * t;
          return innerDeriv * result;
        } else {
          t -= 4.0;
          double result = -6.3131313131313137e-04;
          result = 2.5252525252525255e-03 + result * t;
          result = -3.7878787878787880e-03 + result * t;
          result = 2.5252525252525255e-03 + result * t;
          result = -6.3131313131313137e-04 + result * t;
          return innerDeriv * result;
        }

      case 7:
        if ((t < -7.0) || (t > 7.0)) {
          return 0.0;
        } else if (t < -6.0) {
          t += 7.0;
          double result = -5.7487122884473878e-07;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          return innerDeriv * result;
        } else if (t < -5.0) {
          t += 6.0;
          double result = 7.3008646063281828e-05;
          result = -3.4492273730684327e-06 + result * t;
          result = -8.6230684326710817e-06 + result * t;
          result = -1.1497424576894776e-05 + result * t;
          result = -8.6230684326710817e-06 + result * t;
          result = -3.4492273730684327e-06 + result * t;
          result = -5.7487122884473878e-07 + result * t;
          return innerDeriv * result;
        } else if (t < -4.0) {
          t += 5.0;
          double result = -1.1796357615894040e-03;
          result = 4.3460264900662250e-04 + result * t;
          result = 1.0692604856512142e-03 + result * t;
          result = 1.3796909492273732e-03 + result * t;
          result = 9.6578366445916120e-04 + result * t;
          result = 3.3112582781456954e-04 + result * t;
          result = 3.6791758646063282e-05 + result * t;
          return innerDeriv * result;
        } else if (t < -3.0) {
          t += 4.0;
          double result = 7.6503863134657841e-03;
          result = -6.6432119205298013e-03 + result * t;
          result = -1.4452262693156732e-02 + result * t;
          result = -1.3589955849889625e-02 + result * t;
          result = -1.8280905077262693e-03 + result * t;
          result = 5.7740066225165565e-03 + result * t;
          result = 3.0376195732155996e-03 + result * t;
          return innerDeriv * result;
        } else if (t < -2.0) {
          t += 3.0;
          double result = -2.7219577814569535e-02;
          result = 3.9259105960264898e-02 + result * t;
          result = 6.7087472406181015e-02 + result * t;
          result = 1.5176600441501103e-02 + result * t;
          result = -8.0987858719646796e-02 + result * t;
          result = -8.3774834437086096e-02 + result * t;
          result = -2.0051508462104489e-02 + result * t;
          return innerDeriv * result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = 6.0418391280353201e-02;
          result = -1.2405836092715232e-01 + result * t;
          result = -1.4491066501103753e-01 + result * t;
          result = 1.3172599337748345e-01 + result * t;
          result = 3.5136416942604859e-01 + result * t;
          result = 1.0110720198675496e-01 + result * t;
          result = -9.0510600625459903e-02 + result * t;
          return innerDeriv * result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = -8.8888888888888892e-02;
          result = 2.3845198675496690e-01 + result * t;
          result = 1.4107339955849890e-01 + result * t;
          result = -4.8013245033112584e-01 + result * t;
          result = -4.5722958057395141e-01 + result * t;
          result = 3.6158940397350992e-01 + result * t;
          result = 2.8513612950699041e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 1.0) {
          double result = 8.8888888888888892e-02;
          result = -2.9488134657836645e-01 + result * t;
          result *= t;
          result = 6.9090323767476081e-01 + result * t;
          result *= t;
          result = -7.7004690949227372e-01 + result * t;
          result *= t;
          return innerDeriv * result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = -6.0418391280353201e-02;
          result = 2.3845198675496690e-01 + result * t;
          result = -1.4107339955849890e-01 + result * t;
          result = -4.8013245033112584e-01 + result * t;
          result = 4.5722958057395141e-01 + result * t;
          result = 3.6158940397350992e-01 + result * t;
          result = -2.8513612950699041e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 3.0) {
          t -= 2.0;
          double result = 2.7219577814569535e-02;
          result = -1.2405836092715232e-01 + result * t;
          result = 1.4491066501103753e-01 + result * t;
          result = 1.3172599337748345e-01 + result * t;
          result = -3.5136416942604859e-01 + result * t;
          result = 1.0110720198675496e-01 + result * t;
          result = 9.0510600625459903e-02 + result * t;
          return innerDeriv * result;
        } else if (t < 4.0) {
          t -= 3.0;
          double result = -7.6503863134657841e-03;
          result = 3.9259105960264898e-02 + result * t;
          result = -6.7087472406181015e-02 + result * t;
          result = 1.5176600441501103e-02 + result * t;
          result = 8.0987858719646796e-02 + result * t;
          result = -8.3774834437086096e-02 + result * t;
          result = 2.0051508462104489e-02 + result * t;
          return innerDeriv * result;
        } else if (t < 5.0) {
          t -= 4.0;
          double result = 1.1796357615894040e-03;
          result = -6.6432119205298013e-03 + result * t;
          result = 1.4452262693156732e-02 + result * t;
          result = -1.3589955849889625e-02 + result * t;
          result = 1.8280905077262693e-03 + result * t;
          result = 5.7740066225165565e-03 + result * t;
          result = -3.0376195732155996e-03 + result * t;
          return innerDeriv * result;
        } else if (t < 6.0) {
          t -= 5.0;
          double result = -7.3008646063281828e-05;
          result = 4.3460264900662250e-04 + result * t;
          result = -1.0692604856512142e-03 + result * t;
          result = 1.3796909492273732e-03 + result * t;
          result = -9.6578366445916120e-04 + result * t;
          result = 3.3112582781456954e-04 + result * t;
          result = -3.6791758646063282e-05 + result * t;
          return innerDeriv * result;
        } else {
          t -= 6.0;
          double result = 5.7487122884473878e-07;
          result = -3.4492273730684327e-06 + result * t;
          result = 8.6230684326710817e-06 + result * t;
          result = -1.1497424576894776e-05 + result * t;
          result = 8.6230684326710817e-06 + result * t;
          result = -3.4492273730684327e-06 + result * t;
          result = 5.7487122884473878e-07 + result * t;
          return innerDeriv * result;
        }

      default:
        return 0.0;
    }
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "WeaklyFundamentalSplineBasisDeriv1: evalDx not implemented" << std::endl;
    return -1;
  }

  inline double getIntegral(LT level, IT index) override {
    std::cerr << "WeaklyFundamentalSplineBasisDeriv1: getIntegral not implemented" << std::endl;
    return -1.0;
  }

  /**
   * @return      Spline degree
   */
  inline size_t getDegree() const override { return degree; }

 protected:
  /// degree of the spline
  size_t degree;
};

// default type-def (unsigned int for level and index)
typedef WeaklyFundamentalSplineBasisDeriv1<unsigned int, unsigned int>
    SWeaklyFundamentalSplineBaseDeriv1;

}  // namespace base
}  // namespace sgpp

#endif /* WEAKLY_FUNDAMENTAL_SPLINE_BASE_DERIV1_HPP */
