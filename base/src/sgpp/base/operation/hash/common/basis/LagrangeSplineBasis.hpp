// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LAGRANGE_SPLINE_BASE_HPP
#define LAGRANGE_SPLINE_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <cstddef>
#include <algorithm>

namespace sgpp {
namespace base {

/**
 * Lagrange spline basis.
 */
template <class LT, class IT>
class LagrangeSplineBasis: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  LagrangeSplineBasis() : degree(0) {
  }

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit LagrangeSplineBasis(size_t degree) : degree(degree) {
    if (degree < 1) {
      this->degree = 1;
    } else if (degree % 2 == 0) {
      this->degree = degree - 1;
    }

    if (this->degree > 7) {
      throw std::runtime_error("Unsupported Lagrange spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~LagrangeSplineBasis() override {
  }

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
          double result = -1.0/24.0;
          result *= t;
          result *= t;
          result *= t;
          return result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = 7.0/24.0;
          result = -1.0/8.0 + result * t;
          result = -1.0/8.0 + result * t;
          result = -1.0/24.0 + result * t;
          return result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = -2.0/3.0;
          result = 3.0/4.0 + result * t;
          result = 1.0/2.0 + result * t;
          result *= t;
          return result;
        } else if (t < 1.0) {
          double result = 2.0/3.0;
          result = -5.0/4.0 + result * t;
          result *= t;
          result = 7.0/12.0 + result * t;
          return result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = -7.0/24.0;
          result = 3.0/4.0 + result * t;
          result = -1.0/2.0 + result * t;
          result *= t;
          return result;
        } else {
          t -= 2.0;
          double result = 1.0/24.0;
          result = -1.0/8.0 + result * t;
          result = 1.0/8.0 + result * t;
          result = -1.0/24.0 + result * t;
          return result;
        }

      case 5:
        if ((t < -5.0) || (t > 5.0)) {
          return 0.0;
        } else if (t < -4.0) {
          t += 5.0;
          double result = 1.0/7920.0;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          result *= t;
          return result;
        } else if (t < -3.0) {
          t += 4.0;
          double result = -31.0/7920.0;
          result = 1.0/1584.0 + result * t;
          result = 1.0/792.0 + result * t;
          result = 1.0/792.0 + result * t;
          result = 1.0/1584.0 + result * t;
          result = 1.0/7920.0 + result * t;
          return result;
        } else if (t < -2.0) {
          t += 3.0;
          double result = 103.0/3960.0;
          result = -5.0/264.0 + result * t;
          result = -7.0/198.0 + result * t;
          result = -1.0/33.0 + result * t;
          result = -1.0/99.0 + result * t;
          result *= t;
          return result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = -313.0/3960.0;
          result = 1.0/9.0 + result * t;
          result = 59.0/396.0 + result * t;
          result = 1.0/99.0 + result * t;
          result = -97.0/792.0 + result * t;
          result = -34.0/495.0 + result * t;
          return result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = 2.0/15.0;
          result = -25.0/88.0 + result * t;
          result = -13.0/66.0 + result * t;
          result = 1.0/3.0 + result * t;
          result = 13.0/33.0 + result * t;
          result *= t;
          return result;
        } else if (t < 1.0) {
          double result = -2.0/15.0;
          result = 101.0/264.0 + result * t;
          result *= t;
          result = -83.0/132.0 + result * t;
          result *= t;
          result = 167.0/440.0 + result * t;
          return result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = 313.0/3960.0;
          result = -25.0/88.0 + result * t;
          result = 13.0/66.0 + result * t;
          result = 1.0/3.0 + result * t;
          result = -13.0/33.0 + result * t;
          result *= t;
          return result;
        } else if (t < 3.0) {
          t -= 2.0;
          double result = -103.0/3960.0;
          result = 1.0/9.0 + result * t;
          result = -59.0/396.0 + result * t;
          result = 1.0/99.0 + result * t;
          result = 97.0/792.0 + result * t;
          result = -34.0/495.0 + result * t;
          return result;
        } else if (t < 4.0) {
          t -= 3.0;
          double result = 31.0/7920.0;
          result = -5.0/264.0 + result * t;
          result = 7.0/198.0 + result * t;
          result = -1.0/33.0 + result * t;
          result = 1.0/99.0 + result * t;
          result *= t;
          return result;
        } else {
          t -= 4.0;
          double result = -1.0/7920.0;
          result = 1.0/1584.0 + result * t;
          result = -1.0/792.0 + result * t;
          result = 1.0/792.0 + result * t;
          result = -1.0/1584.0 + result * t;
          result = 1.0/7920.0 + result * t;
          return result;
        }

      case 7:
        if ((t < -7.0) || (t > 7.0)) {
          return 0.0;
        } else if (t < -6.0) {
          t += 7.0;
          double result = -1.0/12176640.0;
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
          double result = 127.0/12176640.0;
          result = -1.0/1739520.0 + result * t;
          result = -1.0/579840.0 + result * t;
          result = -1.0/347904.0 + result * t;
          result = -1.0/347904.0 + result * t;
          result = -1.0/579840.0 + result * t;
          result = -1.0/1739520.0 + result * t;
          result = -1.0/12176640.0 + result * t;
          return result;
        } else if (t < -4.0) {
          t += 5.0;
          double result = -57.0/338240.0;
          result = 7.0/96640.0 + result * t;
          result = 31.0/144960.0 + result * t;
          result = 5.0/14496.0 + result * t;
          result = 7.0/21744.0 + result * t;
          result = 1.0/6040.0 + result * t;
          result = 1.0/27180.0 + result * t;
          result *= t;
          return result;
        } else if (t < -3.0) {
          t += 4.0;
          double result = 1109.0/1014720.0;
          result = -107.0/96640.0 + result * t;
          result = -419.0/144960.0 + result * t;
          result = -197.0/57984.0 + result * t;
          result = -53.0/86976.0 + result * t;
          result = 279.0/96640.0 + result * t;
          result = 1321.0/434880.0 + result * t;
          result = 2003.0/2029440.0 + result * t;
          return result;
        } else if (t < -2.0) {
          t += 3.0;
          double result = -5261.0/1352960.0;
          result = 1897.0/289920.0 + result * t;
          result = 389.0/28992.0 + result * t;
          result = 55.0/14496.0 + result * t;
          result = -587.0/21744.0 + result * t;
          result = -253.0/6040.0 + result * t;
          result = -109.0/5436.0 + result * t;
          result *= t;
          return result;
        } else if (t < -1.0) {
          t += 2.0;
          double result = 35033.0/4058880.0;
          result = -11989.0/579840.0 + result * t;
          result = -3361.0/115968.0 + result * t;
          result = 1273.0/38656.0 + result * t;
          result = 40747.0/347904.0 + result * t;
          result = 9771.0/193280.0 + result * t;
          result = -31489.0/347904.0 + result * t;
          result = -93447.0/1352960.0 + result * t;
          return result;
        } else if (t < 0.0) {
          t += 1.0;
          double result = -4.0/315.0;
          result = 5761.0/144960.0 + result * t;
          result = 409.0/14496.0 + result * t;
          result = -145.0/1208.0 + result * t;
          result = -1657.0/10872.0 + result * t;
          result = 273.0/1510.0 + result * t;
          result = 775.0/2718.0 + result * t;
          result *= t;
          return result;
        } else if (t < 1.0) {
          double result = 4.0/315.0;
          result = -21373.0/434880.0 + result * t;
          result *= t;
          result = 15023.0/86976.0 + result * t;
          result *= t;
          result = -55813.0/144960.0 + result * t;
          result *= t;
          result = 757223.0/3044160.0 + result * t;
          return result;
        } else if (t < 2.0) {
          t -= 1.0;
          double result = -35033.0/4058880.0;
          result = 5761.0/144960.0 + result * t;
          result = -409.0/14496.0 + result * t;
          result = -145.0/1208.0 + result * t;
          result = 1657.0/10872.0 + result * t;
          result = 273.0/1510.0 + result * t;
          result = -775.0/2718.0 + result * t;
          result *= t;
          return result;
        } else if (t < 3.0) {
          t -= 2.0;
          double result = 5261.0/1352960.0;
          result = -11989.0/579840.0 + result * t;
          result = 3361.0/115968.0 + result * t;
          result = 1273.0/38656.0 + result * t;
          result = -40747.0/347904.0 + result * t;
          result = 9771.0/193280.0 + result * t;
          result = 31489.0/347904.0 + result * t;
          result = -93447.0/1352960.0 + result * t;
          return result;
        } else if (t < 4.0) {
          t -= 3.0;
          double result = -1109.0/1014720.0;
          result = 1897.0/289920.0 + result * t;
          result = -389.0/28992.0 + result * t;
          result = 55.0/14496.0 + result * t;
          result = 587.0/21744.0 + result * t;
          result = -253.0/6040.0 + result * t;
          result = 109.0/5436.0 + result * t;
          result *= t;
          return result;
        } else if (t < 5.0) {
          t -= 4.0;
          double result = 57.0/338240.0;
          result = -107.0/96640.0 + result * t;
          result = 419.0/144960.0 + result * t;
          result = -197.0/57984.0 + result * t;
          result = 53.0/86976.0 + result * t;
          result = 279.0/96640.0 + result * t;
          result = -1321.0/434880.0 + result * t;
          result = 2003.0/2029440.0 + result * t;
          return result;
        } else if (t < 6.0) {
          t -= 5.0;
          double result = -127.0/12176640.0;
          result = 7.0/96640.0 + result * t;
          result = -31.0/144960.0 + result * t;
          result = 5.0/14496.0 + result * t;
          result = -7.0/21744.0 + result * t;
          result = 1.0/6040.0 + result * t;
          result = -1.0/27180.0 + result * t;
          result *= t;
          return result;
        } else {
          t -= 6.0;
          double result = 1.0/12176640.0;
          result = -1.0/1739520.0 + result * t;
          result = 1.0/579840.0 + result * t;
          result = -1.0/347904.0 + result * t;
          result = 1.0/347904.0 + result * t;
          result = -1.0/579840.0 + result * t;
          result = 1.0/1739520.0 + result * t;
          result = -1.0/12176640.0 + result * t;
          return result;
        }

      default:
        return 0.0;
    }
  }

  /**
   * @return      Spline degree
   */
  inline size_t getDegree() const {
    return degree;
  }

 protected:
  /// degree of the spline
  size_t degree;
};

// default type-def (unsigned int for level and index)
typedef LagrangeSplineBasis<unsigned int, unsigned int> SLagrangeSplineBase;

}  // namespace base
}  // namespace sgpp

#endif /* LAGRANGE_SPLINE_BASE_HPP */
