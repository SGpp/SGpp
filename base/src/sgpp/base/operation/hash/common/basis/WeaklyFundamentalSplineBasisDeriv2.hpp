// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WEAKLY_FUNDAMENTAL_SPLINE_BASE_DERIV2_HPP
#define WEAKLY_FUNDAMENTAL_SPLINE_BASE_DERIV2_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * weakly fundamental spline basis (2nd derivative).
 */
template <class LT, class IT>
class WeaklyFundamentalSplineBasisDeriv2 : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  WeaklyFundamentalSplineBasisDeriv2() : degree(0) {}

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit WeaklyFundamentalSplineBasisDeriv2(size_t degree) : degree(degree) {
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
  ~WeaklyFundamentalSplineBasisDeriv2() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    const double innerDeriv = hInvDbl * hInvDbl;
    double t = x * hInvDbl - static_cast<double>(i);

    switch (degree) {
      case 1:
        return 0.0;

      case 3:
        if ((t < -3.0) || (t > 3.0)) {
          return 0.0;
        } else if (t < -2.0) {
          t += 3.0;
          double result = -2.5000000000000000e-01;
          result *= t;
          return innerDeriv * result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = 1.7500000000000000e+00;
          result = -2.5000000000000000e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = -4.0;
          result = 1.5000000000000000e+00 + result * t;
          return innerDeriv * result;
        } else if (t < 1.0) {
          double result = 4.0;
          result = -2.5000000000000000e+00 + result * t;
          return innerDeriv * result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = -1.7500000000000000e+00;
          result = 1.5000000000000000e+00 + result * t;
          return innerDeriv * result;
        } else {
          t -= 2.0;
          double result = 2.5000000000000000e-01;
          result = -2.5000000000000000e-01 + result * t;
          return innerDeriv * result;
        }

      case 5:
        if ((t < -5.0) || (t > 5.0)) {
          return 0.0;
        } else if (t < -4.0) {
          t += 5.0;
          double result = 2.5252525252525255e-03;
          result *= t;
          result *= t;
          result *= t;
          return innerDeriv * result;
        } else if (t < -3.0) {
          t += 4.0;
          double result = -7.8282828282828287e-02;
          result = 7.5757575757575760e-03 + result * t;
          result = 7.5757575757575760e-03 + result * t;
          result = 2.5252525252525255e-03 + result * t;
          return innerDeriv * result;
        } else if (t < -2.0) {
          t += 3.0;
          double result = 5.2020202020202022e-01;
          result = -2.2727272727272727e-01 + result * t;
          result = -2.1212121212121213e-01 + result * t;
          result = -6.0606060606060608e-02 + result * t;
          return innerDeriv * result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = -1.5808080808080809e+00;
          result = 1.3333333333333333e+00 + result * t;
          result = 8.9393939393939392e-01 + result * t;
          result = 2.0202020202020204e-02 + result * t;
          return innerDeriv * result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = 2.6666666666666665e+00;
          result = -3.4090909090909092e+00 + result * t;
          result = -1.1818181818181819e+00 + result * t;
          result = 6.6666666666666663e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 1.0) {
          double result = -2.6666666666666665e+00;
          result = 4.5909090909090908e+00 + result * t;
          result *= t;
          result = -1.2575757575757576e+00 + result * t;
          return innerDeriv * result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = 1.5808080808080809e+00;
          result = -3.4090909090909092e+00 + result * t;
          result = 1.1818181818181819e+00 + result * t;
          result = 6.6666666666666663e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 3.0) {
          t -= 2.0;
          double result = -5.2020202020202022e-01;
          result = 1.3333333333333333e+00 + result * t;
          result = -8.9393939393939392e-01 + result * t;
          result = 2.0202020202020204e-02 + result * t;
          return innerDeriv * result;
        } else if (t < 4.0) {
          t -= 3.0;
          double result = 7.8282828282828287e-02;
          result = -2.2727272727272727e-01 + result * t;
          result = 2.1212121212121213e-01 + result * t;
          result = -6.0606060606060608e-02 + result * t;
          return innerDeriv * result;
        } else {
          t -= 4.0;
          double result = -2.5252525252525255e-03;
          result = 7.5757575757575760e-03 + result * t;
          result = -7.5757575757575760e-03 + result * t;
          result = 2.5252525252525255e-03 + result * t;
          return innerDeriv * result;
        }

      case 7:
        if ((t < -7.0) || (t > 7.0)) {
          return 0.0;
        } else if (t < -6.0) {
          t += 7.0;
          double result = -3.4492273730684327e-06;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          return innerDeriv * result;
        } else if (t < -5.0) {
          t += 6.0;
          double result = 4.3805187637969094e-04;
          result = -1.7246136865342163e-05 + result * t;
          result = -3.4492273730684327e-05 + result * t;
          result = -3.4492273730684327e-05 + result * t;
          result = -1.7246136865342163e-05 + result * t;
          result = -3.4492273730684327e-06 + result * t;
          return innerDeriv * result;
        } else if (t < -4.0) {
          t += 5.0;
          double result = -7.0778145695364241e-03;
          result = 2.1730132450331124e-03 + result * t;
          result = 4.2770419426048567e-03 + result * t;
          result = 4.1390728476821195e-03 + result * t;
          result = 1.9315673289183224e-03 + result * t;
          result = 3.3112582781456954e-04 + result * t;
          return innerDeriv * result;
        } else if (t < -3.0) {
          t += 4.0;
          double result = 4.5902317880794699e-02;
          result = -3.3216059602649006e-02 + result * t;
          result = -5.7809050772626928e-02 + result * t;
          result = -4.0769867549668874e-02 + result * t;
          result = -3.6561810154525387e-03 + result * t;
          result = 5.7740066225165565e-03 + result * t;
          return innerDeriv * result;
        } else if (t < -2.0) {
          t += 3.0;
          double result = -1.6331746688741722e-01;
          result = 1.9629552980132450e-01 + result * t;
          result = 2.6834988962472406e-01 + result * t;
          result = 4.5529801324503308e-02 + result * t;
          result = -1.6197571743929359e-01 + result * t;
          result = -8.3774834437086096e-02 + result * t;
          return innerDeriv * result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = 3.6251034768211921e-01;
          result = -6.2029180463576161e-01 + result * t;
          result = -5.7964266004415010e-01 + result * t;
          result = 3.9517798013245031e-01 + result * t;
          result = 7.0272833885209718e-01 + result * t;
          result = 1.0110720198675496e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = -5.3333333333333333e-01;
          result = 1.1922599337748345e+00 + result * t;
          result = 5.6429359823399561e-01 + result * t;
          result = -1.4403973509933774e+00 + result * t;
          result = -9.1445916114790282e-01 + result * t;
          result = 3.6158940397350992e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 1.0) {
          double result = 5.3333333333333333e-01;
          result = -1.4744067328918322e+00 + result * t;
          result *= t;
          result = 2.0727097130242824e+00 + result * t;
          result *= t;
          result = -7.7004690949227372e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = -3.6251034768211921e-01;
          result = 1.1922599337748345e+00 + result * t;
          result = -5.6429359823399561e-01 + result * t;
          result = -1.4403973509933774e+00 + result * t;
          result = 9.1445916114790282e-01 + result * t;
          result = 3.6158940397350992e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 3.0) {
          t -= 2.0;
          double result = 1.6331746688741722e-01;
          result = -6.2029180463576161e-01 + result * t;
          result = 5.7964266004415010e-01 + result * t;
          result = 3.9517798013245031e-01 + result * t;
          result = -7.0272833885209718e-01 + result * t;
          result = 1.0110720198675496e-01 + result * t;
          return innerDeriv * result;
        } else if (t < 4.0) {
          t -= 3.0;
          double result = -4.5902317880794699e-02;
          result = 1.9629552980132450e-01 + result * t;
          result = -2.6834988962472406e-01 + result * t;
          result = 4.5529801324503308e-02 + result * t;
          result = 1.6197571743929359e-01 + result * t;
          result = -8.3774834437086096e-02 + result * t;
          return innerDeriv * result;
        } else if (t < 5.0) {
          t -= 4.0;
          double result = 7.0778145695364241e-03;
          result = -3.3216059602649006e-02 + result * t;
          result = 5.7809050772626928e-02 + result * t;
          result = -4.0769867549668874e-02 + result * t;
          result = 3.6561810154525387e-03 + result * t;
          result = 5.7740066225165565e-03 + result * t;
          return innerDeriv * result;
        } else if (t < 6.0) {
          t -= 5.0;
          double result = -4.3805187637969094e-04;
          result = 2.1730132450331124e-03 + result * t;
          result = -4.2770419426048567e-03 + result * t;
          result = 4.1390728476821195e-03 + result * t;
          result = -1.9315673289183224e-03 + result * t;
          result = 3.3112582781456954e-04 + result * t;
          return innerDeriv * result;
        } else {
          t -= 6.0;
          double result = 3.4492273730684327e-06;
          result = -1.7246136865342163e-05 + result * t;
          result = 3.4492273730684327e-05 + result * t;
          result = -3.4492273730684327e-05 + result * t;
          result = 1.7246136865342163e-05 + result * t;
          result = -3.4492273730684327e-06 + result * t;
          return innerDeriv * result;
        }

      default:
        return 0.0;
    }
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "WeaklyFundamentalSplineBasisDeriv2: evalDx not implemented" << std::endl;
    return -1;
  }

  inline double getIntegral(LT level, IT index) override {
    std::cerr << "WeaklyFundamentalSplineBasisDeriv2: getIntegral not implemented" << std::endl;
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
typedef WeaklyFundamentalSplineBasisDeriv2<unsigned int, unsigned int>
    SWeaklyFundamentalSplineBaseDeriv2;

}  // namespace base
}  // namespace sgpp

#endif /* WEAKLY_FUNDAMENTAL_SPLINE_BASE_DERIV2_HPP */
