// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WEAKLY_FUNDAMENTAL_SPLINE_BASE_HPP
#define WEAKLY_FUNDAMENTAL_SPLINE_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * weakly fundamental spline basis.
 */
template <class LT, class IT>
class WeaklyFundamentalSplineBasis : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  WeaklyFundamentalSplineBasis() : degree(0) {}

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit WeaklyFundamentalSplineBasis(size_t degree) : degree(degree) {
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
  ~WeaklyFundamentalSplineBasis() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    double t = x * static_cast<double>(hInv) - static_cast<double>(i);

    switch (degree) {
      case 1:
        return std::max(1.0 - std::abs(t), 0.0);

      case 3:
        if ((t < -3.0) || (t > 3.0)) {
          return 0.0;
        } else if (t < -2.0) {
          t += 3.0;
          double result = -4.1666666666666664e-02;
          result *= t;
          result *= t;
          result *= t;
          return result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = 2.9166666666666669e-01;
          result = -1.2500000000000000e-01 + result * t;
          result = -1.2500000000000000e-01 + result * t;
          result = -4.1666666666666664e-02 + result * t;
          return result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = -6.6666666666666663e-01;
          result = 7.5000000000000000e-01 + result * t;
          result = 5.0000000000000000e-01 + result * t;
          result *= t;
          return result;
        } else if (t < 1.0) {
          double result = 6.6666666666666663e-01;
          result = -1.2500000000000000e+00 + result * t;
          result *= t;
          result = 5.8333333333333337e-01 + result * t;
          return result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = -2.9166666666666669e-01;
          result = 7.5000000000000000e-01 + result * t;
          result = -5.0000000000000000e-01 + result * t;
          result *= t;
          return result;
        } else {
          t -= 2.0;
          double result = 4.1666666666666664e-02;
          result = -1.2500000000000000e-01 + result * t;
          result = 1.2500000000000000e-01 + result * t;
          result = -4.1666666666666664e-02 + result * t;
          return result;
        }

      case 5:
        if ((t < -5.0) || (t > 5.0)) {
          return 0.0;
        } else if (t < -4.0) {
          t += 5.0;
          double result = 1.2626262626262626e-04;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          return result;
        } else if (t < -3.0) {
          t += 4.0;
          double result = -3.9141414141414142e-03;
          result = 6.3131313131313137e-04 + result * t;
          result = 1.2626262626262627e-03 + result * t;
          result = 1.2626262626262627e-03 + result * t;
          result = 6.3131313131313137e-04 + result * t;
          result = 1.2626262626262626e-04 + result * t;
          return result;
        } else if (t < -2.0) {
          t += 3.0;
          double result = 2.6010101010101011e-02;
          result = -1.8939393939393940e-02 + result * t;
          result = -3.5353535353535352e-02 + result * t;
          result = -3.0303030303030304e-02 + result * t;
          result = -1.0101010101010102e-02 + result * t;
          result *= t;
          return result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = -7.9040404040404036e-02;
          result = 1.1111111111111110e-01 + result * t;
          result = 1.4898989898989898e-01 + result * t;
          result = 1.0101010101010102e-02 + result * t;
          result = -1.2247474747474747e-01 + result * t;
          result = -6.8686868686868685e-02 + result * t;
          return result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = 1.3333333333333333e-01;
          result = -2.8409090909090912e-01 + result * t;
          result = -1.9696969696969696e-01 + result * t;
          result = 3.3333333333333331e-01 + result * t;
          result = 3.9393939393939392e-01 + result * t;
          result *= t;
          return result;
        } else if (t < 1.0) {
          double result = -1.3333333333333333e-01;
          result = 3.8257575757575757e-01 + result * t;
          result *= t;
          result = -6.2878787878787878e-01 + result * t;
          result *= t;
          result = 3.7954545454545452e-01 + result * t;
          return result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = 7.9040404040404036e-02;
          result = -2.8409090909090912e-01 + result * t;
          result = 1.9696969696969696e-01 + result * t;
          result = 3.3333333333333331e-01 + result * t;
          result = -3.9393939393939392e-01 + result * t;
          result *= t;
          return result;
        } else if (t < 3.0) {
          t -= 2.0;
          double result = -2.6010101010101011e-02;
          result = 1.1111111111111110e-01 + result * t;
          result = -1.4898989898989898e-01 + result * t;
          result = 1.0101010101010102e-02 + result * t;
          result = 1.2247474747474747e-01 + result * t;
          result = -6.8686868686868685e-02 + result * t;
          return result;
        } else if (t < 4.0) {
          t -= 3.0;
          double result = 3.9141414141414142e-03;
          result = -1.8939393939393940e-02 + result * t;
          result = 3.5353535353535352e-02 + result * t;
          result = -3.0303030303030304e-02 + result * t;
          result = 1.0101010101010102e-02 + result * t;
          result *= t;
          return result;
        } else {
          t -= 4.0;
          double result = -1.2626262626262626e-04;
          result = 6.3131313131313137e-04 + result * t;
          result = -1.2626262626262627e-03 + result * t;
          result = 1.2626262626262627e-03 + result * t;
          result = -6.3131313131313137e-04 + result * t;
          result = 1.2626262626262626e-04 + result * t;
          return result;
        }

      case 7:
        if ((t < -7.0) || (t > 7.0)) {
          return 0.0;
        } else if (t < -6.0) {
          t += 7.0;
          double result = -8.2124461263534111e-08;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          return result;
        } else if (t < -5.0) {
          t += 6.0;
          double result = 1.0429806580468832e-05;
          result = -5.7487122884473878e-07 + result * t;
          result = -1.7246136865342163e-06 + result * t;
          result = -2.8743561442236939e-06 + result * t;
          result = -2.8743561442236939e-06 + result * t;
          result = -1.7246136865342163e-06 + result * t;
          result = -5.7487122884473878e-07 + result * t;
          result = -8.2124461263534111e-08 + result * t;
          return result;
        } else if (t < -4.0) {
          t += 5.0;
          double result = -1.6851939451277201e-04;
          result = 7.2433774834437093e-05 + result * t;
          result = 2.1385209713024284e-04 + result * t;
          result = 3.4492273730684329e-04 + result * t;
          result = 3.2192788815305372e-04 + result * t;
          result = 1.6556291390728477e-04 + result * t;
          result = 3.6791758646063282e-05 + result * t;
          result *= t;
          return result;
        } else if (t < -3.0) {
          t += 4.0;
          double result = 1.0929123304951119e-03;
          result = -1.1072019867549669e-03 + result * t;
          result = -2.8904525386313465e-03 + result * t;
          result = -3.3974889624724062e-03 + result * t;
          result = -6.0936350257542315e-04 + result * t;
          result = 2.8870033112582782e-03 + result * t;
          result = 3.0376195732155996e-03 + result * t;
          result = 9.8697177546515293e-04 + result * t;
          return result;
        } else if (t < -2.0) {
          t += 3.0;
          double result = -3.8885111163670764e-03;
          result = 6.5431843267108169e-03 + result * t;
          result = 1.3417494481236202e-02 + result * t;
          result = 3.7941501103752758e-03 + result * t;
          result = -2.6995952906548933e-02 + result * t;
          result = -4.1887417218543048e-02 + result * t;
          result = -2.0051508462104489e-02 + result * t;
          result *= t;
          return result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = 8.6311987543361713e-03;
          result = -2.0676393487858720e-02 + result * t;
          result = -2.8982133002207505e-02 + result * t;
          result = 3.2931498344370862e-02 + result * t;
          result = 1.1712138980868285e-01 + result * t;
          result = 5.0553600993377482e-02 + result * t;
          result = -9.0510600625459903e-02 + result * t;
          result = -6.9068560785241248e-02 + result * t;
          return result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = -1.2698412698412698e-02;
          result = 3.9741997792494481e-02 + result * t;
          result = 2.8214679911699778e-02 + result * t;
          result = -1.2003311258278146e-01 + result * t;
          result = -1.5240986019131714e-01 + result * t;
          result = 1.8079470198675496e-01 + result * t;
          result = 2.8513612950699041e-01 + result * t;
          result *= t;
          return result;
        } else if (t < 1.0) {
          double result = 1.2698412698412698e-02;
          result = -4.9146891096394404e-02 + result * t;
          result *= t;
          result = 1.7272580941869020e-01 + result * t;
          result *= t;
          result = -3.8502345474613686e-01 + result * t;
          result *= t;
          result = 2.4874612372542837e-01 + result * t;
          return result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = -8.6311987543361713e-03;
          result = 3.9741997792494481e-02 + result * t;
          result = -2.8214679911699778e-02 + result * t;
          result = -1.2003311258278146e-01 + result * t;
          result = 1.5240986019131714e-01 + result * t;
          result = 1.8079470198675496e-01 + result * t;
          result = -2.8513612950699041e-01 + result * t;
          result *= t;
          return result;
        } else if (t < 3.0) {
          t -= 2.0;
          double result = 3.8885111163670764e-03;
          result = -2.0676393487858720e-02 + result * t;
          result = 2.8982133002207505e-02 + result * t;
          result = 3.2931498344370862e-02 + result * t;
          result = -1.1712138980868285e-01 + result * t;
          result = 5.0553600993377482e-02 + result * t;
          result = 9.0510600625459903e-02 + result * t;
          result = -6.9068560785241248e-02 + result * t;
          return result;
        } else if (t < 4.0) {
          t -= 3.0;
          double result = -1.0929123304951119e-03;
          result = 6.5431843267108169e-03 + result * t;
          result = -1.3417494481236202e-02 + result * t;
          result = 3.7941501103752758e-03 + result * t;
          result = 2.6995952906548933e-02 + result * t;
          result = -4.1887417218543048e-02 + result * t;
          result = 2.0051508462104489e-02 + result * t;
          result *= t;
          return result;
        } else if (t < 5.0) {
          t -= 4.0;
          double result = 1.6851939451277201e-04;
          result = -1.1072019867549669e-03 + result * t;
          result = 2.8904525386313465e-03 + result * t;
          result = -3.3974889624724062e-03 + result * t;
          result = 6.0936350257542315e-04 + result * t;
          result = 2.8870033112582782e-03 + result * t;
          result = -3.0376195732155996e-03 + result * t;
          result = 9.8697177546515293e-04 + result * t;
          return result;
        } else if (t < 6.0) {
          t -= 5.0;
          double result = -1.0429806580468832e-05;
          result = 7.2433774834437093e-05 + result * t;
          result = -2.1385209713024284e-04 + result * t;
          result = 3.4492273730684329e-04 + result * t;
          result = -3.2192788815305372e-04 + result * t;
          result = 1.6556291390728477e-04 + result * t;
          result = -3.6791758646063282e-05 + result * t;
          result *= t;
          return result;
        } else {
          t -= 6.0;
          double result = 8.2124461263534111e-08;
          result = -5.7487122884473878e-07 + result * t;
          result = 1.7246136865342163e-06 + result * t;
          result = -2.8743561442236939e-06 + result * t;
          result = 2.8743561442236939e-06 + result * t;
          result = -1.7246136865342163e-06 + result * t;
          result = 5.7487122884473878e-07 + result * t;
          result = -8.2124461263534111e-08 + result * t;
          return result;
        }

      default:
        return 0.0;
    }
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "WeaklyFundamentalSplineBasis: evalDx not implemented" << std::endl;
    return -1;
  }

  inline double getIntegral(LT level, IT index) override {
    std::cerr << "WeaklyFundamentalSplineBasis: getIntegral not implemented" << std::endl;
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
typedef WeaklyFundamentalSplineBasis<unsigned int, unsigned int> SWeaklyFundamentalSplineBase;

}  // namespace base
}  // namespace sgpp

#endif /* WEAKLY_FUNDAMENTAL_SPLINE_BASE_HPP */
