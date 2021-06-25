// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NATURAL_BSPLINE_BASE_HPP
#define NATURAL_BSPLINE_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * B-spline basis with natural boundary conditions.
 */
template <class LT, class IT>
class NaturalBsplineBasis : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NaturalBsplineBasis() : bsplineBasis(BsplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NaturalBsplineBasis(size_t degree) : bsplineBasis(BsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NaturalBsplineBasis() override {}

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
        return std::max(1.0 - std::abs(t), 0.0);

      case 3:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return 1.0 - x;
          } else {
            // l = 0, i = 1
            return x;
          }
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return bsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 1) {
            // l = 1, i = 1
            if ((t < -1.0) || (t > 1.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -1.0 / 6.0;
              result *= t;
              result = 1.0 / 2.0 + result * t;
              result = 1.0 / 2.0 + result * t;
              return result;
            } else {
              double result = 1.0 / 6.0;
              result = -1.0 / 2.0 + result * t;
              result *= t;
              result = 5.0 / 6.0 + result * t;
              return result;
            }
          } else {
            // l >= 2, i = 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -1.0 / 4.0;
              result *= t;
              result = 1.0 / 2.0 + result * t;
              result = 1.0 / 2.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 5.0 / 12.0;
              result = -3.0 / 4.0 + result * t;
              result = -1.0 / 4.0 + result * t;
              result = 3.0 / 4.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0 / 6.0;
              result = 1.0 / 2.0 + result * t;
              result = -1.0 / 2.0 + result * t;
              result = 1.0 / 6.0 + result * t;
              return result;
            }
          }
        }

      case 5:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return 1.0 - x;
          } else {
            // l = 0, i = 1
            return x;
          }
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return bsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 1) {
            // l = 1, i = 1
            if ((t < -1.0) || (t > 1.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 1.0 / 30.0;
              result = -1.0 / 12.0 + result * t;
              result *= t;
              result *= t;
              result = 1.0 / 6.0 + result * t;
              result = 11.0 / 15.0 + result * t;
              return result;
            } else {
              double result = -1.0 / 30.0;
              result = 1.0 / 12.0 + result * t;
              result *= t;
              result = -1.0 / 6.0 + result * t;
              result *= t;
              result = 17.0 / 20.0 + result * t;
              return result;
            }
          } else {
            // l >= 2, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 5.0 / 72.0;
              result = -1.0 / 6.0 + result * t;
              result *= t;
              result *= t;
              result *= t;
              result = 4.0 / 5.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -23.0 / 360.0;
              result = 13.0 / 72.0 + result * t;
              result = 1.0 / 36.0 + result * t;
              result = -11.0 / 36.0 + result * t;
              result = -23.0 / 72.0 + result * t;
              result = 253.0 / 360.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 13.0 / 360.0;
              result = -5.0 / 36.0 + result * t;
              result = 1.0 / 9.0 + result * t;
              result = 2.0 / 9.0 + result * t;
              result = -4.0 / 9.0 + result * t;
              result = 2.0 / 9.0 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.0 / 120.0;
              result = 1.0 / 24.0 + result * t;
              result = -1.0 / 12.0 + result * t;
              result = 1.0 / 12.0 + result * t;
              result = -1.0 / 24.0 + result * t;
              result = 1.0 / 120.0 + result * t;
              return result;
            }
          }
        }

      case 7:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return 1.0 - x;
          } else {
            // l = 0, i = 1
            return x;
          }
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 4, 3 < i < 2^l - 3
          return bsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 1) {
            // l = 1, i = 1
            if ((t < -1.0) || (t > 1.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 1.0 / 840.0;
              result = -1.0 / 192.0 + result * t;
              result = 1.0 / 160.0 + result * t;
              result *= t;
              result *= t;
              result *= t;
              result = -1.0 / 120.0 + result * t;
              result = 29.0 / 28.0 + result * t;
              return result;
            } else {
              double result = -1.0 / 840.0;
              result = 1.0 / 320.0 + result * t;
              result *= t;
              result = -1.0 / 192.0 + result * t;
              result *= t;
              result = 3.0 / 320.0 + result * t;
              result *= t;
              result = 6919.0 / 6720.0 + result * t;
              return result;
            }
          } else if (l == 2) {
            // l = 2, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 5983.0 / 1632960.0;
              result = -205.0 / 7776.0 + result * t;
              result = 59.0 / 1080.0 + result * t;
              result *= t;
              result *= t;
              result *= t;
              result = -209.0 / 135.0 + result * t;
              result = 502.0 / 189.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 3103.0 / 1632960.0;
              result = -167.0 / 233280.0 + result * t;
              result = -2069.0 / 77760.0 + result * t;
              result = 277.0 / 46656.0 + result * t;
              result = 6871.0 / 46656.0 + result * t;
              result = 17713.0 / 77760.0 + result * t;
              result = -328349.0 / 233280.0 + result * t;
              result = 1861357.0 / 1632960.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -5753.0 / 1632960.0;
              result = 367.0 / 29160.0 + result * t;
              result = 35.0 / 3888.0 + result * t;
              result = -52.0 / 729.0 + result * t;
              result = -125.0 / 2916.0 + result * t;
              result = 569.0 / 1215.0 + result * t;
              result = -445.0 / 729.0 + result * t;
              result = 2243.0 / 25515.0 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = 463.0 / 181440.0;
              result = -313.0 / 25920.0 + result * t;
              result = 91.0 / 8640.0 + result * t;
              result = 203.0 / 5184.0 + result * t;
              result = -569.0 / 5184.0 + result * t;
              result = 1007.0 / 8640.0 + result * t;
              result = 211.0 / 25920.0 + result * t;
              result = -27277.0 / 181440.0 + result * t;
              return result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 179.0 / 51840.0;
              result = -19.0 / 864.0 + result * t;
              result = 1.0 / 24.0 + result * t;
              result *= t;
              result *= t;
              result *= t;
              result = -7.0 / 6.0 + result * t;
              result = 95.0 / 42.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 173.0 / 362880.0;
              result = 113.0 / 51840.0 + result * t;
              result = -307.0 / 17280.0 + result * t;
              result = -7.0 / 10368.0 + result * t;
              result = 1013.0 / 10368.0 + result * t;
              result = 2753.0 / 17280.0 + result * t;
              result = -55267.0 / 51840.0 + result * t;
              result = 405833.0 / 362880.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -547.0 / 362880.0;
              result = 143.0 / 25920.0 + result * t;
              result = 23.0 / 4320.0 + result * t;
              result = -13.0 / 324.0 + result * t;
              result = -29.0 / 1296.0 + result * t;
              result = 677.0 / 2160.0 + result * t;
              result = -3431.0 / 6480.0 + result * t;
              result = 13313.0 / 45360.0 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 37.0 / 40320.0;
              result = -29.0 / 5760.0 + result * t;
              result = 13.0 / 1920.0 + result * t;
              result = 19.0 / 1152.0 + result * t;
              result = -83.0 / 1152.0 + result * t;
              result = 211.0 / 1920.0 + result * t;
              result = -467.0 / 5760.0 + result * t;
              result = 979.0 / 40320.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0 / 5040.0;
              result = 1.0 / 720.0 + result * t;
              result = -1.0 / 240.0 + result * t;
              result = 1.0 / 144.0 + result * t;
              result = -1.0 / 144.0 + result * t;
              result = 1.0 / 240.0 + result * t;
              result = -1.0 / 720.0 + result * t;
              result = 1.0 / 5040.0 + result * t;
              return result;
            }
          } else {
            // l >= 3, i = 3
            if ((t < -3.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < -2.0) {
              t += 3.0;
              double result = -161.0 / 51840.0;
              result = 5.0 / 864.0 + result * t;
              result = 1.0 / 120.0 + result * t;
              result *= t;
              result *= t;
              result *= t;
              result = -1.0 / 15.0 + result * t;
              result = -2.0 / 21.0 + result * t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = 343.0 / 51840.0;
              result = -827.0 / 51840.0 + result * t;
              result = -383.0 / 17280.0 + result * t;
              result = 205.0 / 10368.0 + result * t;
              result = 937.0 / 10368.0 + result * t;
              result = 1813.0 / 17280.0 + result * t;
              result = -623.0 / 51840.0 + result * t;
              result = -10951.0 / 72576.0 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -97.0 / 10368.0;
              result = 787.0 / 25920.0 + result * t;
              result = 91.0 / 4320.0 + result * t;
              result = -8.0 / 81.0 + result * t;
              result = -181.0 / 1296.0 + result * t;
              result = 373.0 / 2160.0 + result * t;
              result = 2513.0 / 6480.0 + result * t;
              result = 937.0 / 45360.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 3061.0 / 362880.0;
              result = -607.0 / 17280.0 + result * t;
              result = 13.0 / 1920.0 + result * t;
              result = 155.0 / 1152.0 + result * t;
              result = -17.0 / 384.0 + result * t;
              result = -709.0 / 1920.0 + result * t;
              result = 271.0 / 1920.0 + result * t;
              result = 3103.0 / 8064.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -85.0 / 18144.0;
              result = 31.0 / 1296.0 + result * t;
              result = -29.0 / 1080.0 + result * t;
              result = -41.0 / 648.0 + result * t;
              result = 25.0 / 162.0 + result * t;
              result = 5.0 / 216.0 + result * t;
              result = -251.0 / 810.0 + result * t;
              result = 5149.0 / 22680.0 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = 19.0 / 12960.0;
              result = -23.0 / 2592.0 + result * t;
              result = 79.0 / 4320.0 + result * t;
              result = -7.0 / 2592.0 + result * t;
              result = -137.0 / 2592.0 + result * t;
              result = 85.0 / 864.0 + result * t;
              result = -1001.0 / 12960.0 + result * t;
              result = 2153.0 / 90720.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0 / 5040.0;
              result = 1.0 / 720.0 + result * t;
              result = -1.0 / 240.0 + result * t;
              result = 1.0 / 144.0 + result * t;
              result = -1.0 / 144.0 + result * t;
              result = 1.0 / 240.0 + result * t;
              result = -1.0 / 720.0 + result * t;
              result = 1.0 / 5040.0 + result * t;
              return result;
            }
          }
        }

      default:
        return 0.0;
    }
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "NaturalBsplineBasis: evalDx not implemented" << std::endl;
    return -1;
  }

  inline double getIntegral(LT level, IT index) override {
    std::cerr << "NaturalBsplineBasis: getIntegral not implemented" << std::endl;
    return -1.0;
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override { return bsplineBasis.getDegree(); }

 protected:
  /// B-spline basis for B-spline evaluation
  BsplineBasis<LT, IT> bsplineBasis;
};

// default type-def (unsigned int for level and index)
typedef NaturalBsplineBasis<unsigned int, unsigned int> SNaturalBsplineBase;

}  // namespace base
}  // namespace sgpp

#endif /* NATURAL_BSPLINE_BASE_HPP */
