// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NAK_BSPLINE_MODIFIED_BASE_HPP
#define NAK_BSPLINE_MODIFIED_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>

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
class NakBsplineModifiedBasis: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakBsplineModifiedBasis() : nakBsplineBasis(NakBsplineBasis<LT, IT>()) {
  }

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineModifiedBasis(size_t degree) :
        nakBsplineBasis(NakBsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakBsplineModifiedBasis() override {
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

    switch (getDegree()) {
      case 1:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return std::max(1.0 - std::abs(t), 0.0);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          // l >= 3, i = 1
          return std::max(1.0 - t, 0.0);
        }

      case 3:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return nakBsplineBasis.eval(l, i, x);
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
              double result = 1.0/40.0;
              result *= t;
              result = -3.0/5.0 + result * t;
              result = 6.0/5.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0/40.0;
              result = 3.0/20.0 + result * t;
              result = -3.0/10.0 + result * t;
              result = 1.0/5.0 + result * t;
              return result;
            }
          } else {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 1.0/24.0;
              result *= t;
              result = -3.0/4.0 + result * t;
              result = 5.0/4.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0/12.0;
              result = 1.0/4.0 + result * t;
              result = -1.0/4.0 + result * t;
              result = 1.0/12.0 + result * t;
              return result;
            }
          }
        }

      case 5:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return nakBsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -3.0/70.0;
            result = 19.0/210.0 + result * t;
            result = 37.0/70.0 + result * t;
            result = -331.0/210.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 37.0/22680.0;
              result = -1.0/54.0 + result * t;
              result = 4.0/63.0 + result * t;
              result *= t;
              result = -8.0/21.0 + result * t;
              result = 8.0/15.0 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.0/840.0;
              result = 1.0/168.0 + result * t;
              result = -1.0/84.0 + result * t;
              result = 1.0/84.0 + result * t;
              result = -1.0/168.0 + result * t;
              result = 1.0/840.0 + result * t;
              return result;
            }
          }
        }

      case 7:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return nakBsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -3.0/70.0;
            result = 19.0/210.0 + result * t;
            result = 37.0/70.0 + result * t;
            result = -331.0/210.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else if (l == 3) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 31.0/4942080.0;
              result = -7.0/41184.0 + result * t;
              result = 49.0/25740.0 + result * t;
              result = -14.0/1287.0 + result * t;
              result = 112.0/3861.0 + result * t;
              result *= t;
              result = -3584.0/19305.0 + result * t;
              result = 2048.0/6435.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0/4942080.0;
              result = 7.0/1235520.0 + result * t;
              result = -7.0/102960.0 + result * t;
              result = 7.0/15444.0 + result * t;
              result = -7.0/3861.0 + result * t;
              result = 28.0/6435.0 + result * t;
              result = -112.0/19305.0 + result * t;
              result = 64.0/19305.0 + result * t;
              return result;
            }
          } else {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 2101.0/154828800.0;
              result = -41.0/122880.0 + result * t;
              result = 61.0/18432.0 + result * t;
              result = -25.0/1536.0 + result * t;
              result = 125.0/3456.0 + result * t;
              result *= t;
              result = -125.0/864.0 + result * t;
              result = 125.0/672.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0/151200.0;
              result = 1.0/21600.0 + result * t;
              result = -1.0/7200.0 + result * t;
              result = 1.0/4320.0 + result * t;
              result = -1.0/4320.0 + result * t;
              result = 1.0/7200.0 + result * t;
              result = -1.0/21600.0 + result * t;
              result = 1.0/151200.0 + result * t;
              return result;
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
  inline size_t getDegree() const override {
    return nakBsplineBasis.getDegree();
  }

 protected:
  /// B-spline basis for B-spline evaluation
  NakBsplineBasis<LT, IT> nakBsplineBasis;
};

// default type-def (unsigned int for level and index)
typedef NakBsplineModifiedBasis<unsigned int, unsigned int> SNakBsplineModifiedBase;

}  // namespace base
}  // namespace sgpp

#endif /* NAK_BSPLINE_MODIFIED_BASE_HPP */
