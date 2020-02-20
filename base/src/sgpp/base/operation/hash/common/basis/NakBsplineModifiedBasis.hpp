// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>

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
class NakBsplineModifiedBasis : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakBsplineModifiedBasis() : nakBsplineBasis(NakBsplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineModifiedBasis(size_t degree)
      : nakBsplineBasis(NakBsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakBsplineModifiedBasis() override {}

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

          // l >= 2, i = 1
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
              double result = 1.0 / 40.0;
              result *= t;
              result = -3.0 / 5.0 + result * t;
              result = 6.0 / 5.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0 / 40.0;
              result = 3.0 / 20.0 + result * t;
              result = -3.0 / 10.0 + result * t;
              result = 1.0 / 5.0 + result * t;
              return result;
            }
          } else {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 1.0 / 24.0;
              result *= t;
              result = -3.0 / 4.0 + result * t;
              result = 5.0 / 4.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0 / 12.0;
              result = 1.0 / 4.0 + result * t;
              result = -1.0 / 4.0 + result * t;
              result = 1.0 / 12.0 + result * t;
              return result;
            }
          }
        }

      case 5:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 3, 3 < i < 2^l - 3
          return nakBsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1 (Lagrange)
            double result = -1.0 / 46.0;
            result = 2.0 / 23.0 + result * t;
            result = 9.0 / 23.0 + result * t;
            result = -67.0 / 46.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else if ((l == 3) && (i == 3)) {
            // l = 3, i=3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 2.7557319223985889e-04;
              result = 6.8783068783068784e-03 + result * t;
              result = -5.6084656084656084e-02 + result * t;
              result *= t;
              result = 3.4973544973544973e-01 + result * t;
              result = 3.8603174603174606e-01 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -7.2552910052910051e-03;
              result = 1.1011904761904763e-02 + result * t;
              result = 5.1256613756613757e-02 + result * t;
              result = -5.8928571428571427e-02 + result * t;
              result = -3.1008597883597883e-01 + result * t;
              result = 5.4505952380952383e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 5.7473544973544975e-03;
              result = -2.5264550264550264e-02 + result * t;
              result = 2.2751322751322751e-02 + result * t;
              result = 8.8359788359788360e-02 + result * t;
              result = -2.6640211640211642e-01 + result * t;
              result = 2.3105820105820105e-01 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -2.3148148148148149e-04;
              result = 3.4722222222222220e-03 + result * t;
              result = -2.0833333333333332e-02 + result * t;
              result = 6.2500000000000000e-02 + result * t;
              result = -9.3750000000000000e-02 + result * t;
              result = 5.6250000000000001e-02 + result * t;
              return result;
            }
          } else {
            if (i == 1) {
              // l >= 3, i = 1
              if ((t < -1.0) || (t > 3.0)) {
                return 0.0;
              } else if (t < 2.0) {
                t += 1.0;
                double result = 1.6313932980599648e-03;
                result = -1.8518518518518517e-02 + result * t;
                result = 6.3492063492063489e-02 + result * t;
                result *= t;
                result = -3.8095238095238093e-01 + result * t;
                result = 5.3333333333333333e-01 + result * t;
                return result;
              } else {
                t -= 2.0;
                double result = -1.1904761904761906e-03;
                result = 5.9523809523809521e-03 + result * t;
                result = -1.1904761904761904e-02 + result * t;
                result = 1.1904761904761904e-02 + result * t;
                result = -5.9523809523809521e-03 + result * t;
                result = 1.1904761904761906e-03 + result * t;
                return result;
              }
            } else {
              // l >= 4, i = 3
              if ((t < -3.0) || (t > 3.0)) {
                return 0.0;
              } else if (t < 0.0) {
                t += 3.0;
                double result = 7.9365079365079365e-04;
                result = 4.7619047619047623e-03 + result * t;
                result = -5.7142857142857141e-02 + result * t;
                result *= t;
                result = 3.4285714285714286e-01 + result * t;
                result = 3.7714285714285717e-01 + result * t;
                return result;
              } else if (t < 1.0) {
                double result = -1.2539682539682540e-02;
                result = 1.6666666666666666e-02 + result * t;
                result = 7.1428571428571425e-02 + result * t;
                result = -4.2857142857142858e-02 + result * t;
                result = -3.6428571428571427e-01 + result * t;
                result = 4.4142857142857145e-01 + result * t;
                return result;
              } else if (t < 2.0) {
                t -= 1.0;
                double result = 1.3174603174603174e-02;
                result = -4.6031746031746035e-02 + result * t;
                result = 1.2698412698412698e-02 + result * t;
                result = 1.4603174603174604e-01 + result * t;
                result = -2.3174603174603176e-01 + result * t;
                result = 1.0984126984126984e-01 + result * t;
                return result;
              } else {
                t -= 2.0;
                double result = -3.9682539682539680e-03;
                result = 1.9841269841269840e-02 + result * t;
                result = -3.9682539682539680e-02 + result * t;
                result = 3.9682539682539680e-02 + result * t;
                result = -1.9841269841269840e-02 + result * t;
                result = 3.9682539682539680e-03 + result * t;
                return result;
              }
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
  inline size_t getDegree() const override { return nakBsplineBasis.getDegree(); }

 protected:
  /// B-spline basis for B-spline evaluation
  NakBsplineBasis<LT, IT> nakBsplineBasis;
};

// default type-def (unsigned int for level and index)
typedef NakBsplineModifiedBasis<unsigned int, unsigned int> SNakBsplineModifiedBase;

}  // namespace base
}  // namespace sgpp
