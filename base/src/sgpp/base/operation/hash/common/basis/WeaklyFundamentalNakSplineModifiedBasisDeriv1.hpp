// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WEAKLY_FUNDAMENTAL_NAK_SPLINE_MODIFIED_BASE_DERIV1_HPP
#define WEAKLY_FUNDAMENTAL_NAK_SPLINE_MODIFIED_BASE_DERIV1_HPP

#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineBasisDeriv1.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * weakly fundamental not-a-knot spline basis.
 */
template <class LT, class IT>
class WeaklyFundamentalNakSplineModifiedBasisDeriv1 : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  WeaklyFundamentalNakSplineModifiedBasisDeriv1()
      : weaklyFundamentalNakSplineBasisDeriv1(WeaklyFundamentalNakSplineBasisDeriv1<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit WeaklyFundamentalNakSplineModifiedBasisDeriv1(size_t degree)
      : weaklyFundamentalNakSplineBasisDeriv1(
            WeaklyFundamentalNakSplineBasisDeriv1<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~WeaklyFundamentalNakSplineModifiedBasisDeriv1() override {}

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
        return weaklyFundamentalNakSplineBasisDeriv1.eval(l, i, x);

      case 3:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return weaklyFundamentalNakSplineBasisDeriv1.eval(l, i, x);
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
              double result = 1.2980769230769232e-01;
              result *= t;
              result = -8.0769230769230771e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -1.7307692307692307e-01;
              result = 5.1923076923076927e-01 + result * t;
              result = -2.8846153846153844e-01 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 1.6071428571428573e-01;
              result *= t;
              result = -8.5714285714285710e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -3.7500000000000000e-01;
              result = 6.4285714285714290e-01 + result * t;
              result = -2.1428571428571427e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = 5.3571428571428568e-02;
              result = -1.0714285714285714e-01 + result * t;
              result = 5.3571428571428568e-02 + result * t;
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
          return weaklyFundamentalNakSplineBasisDeriv1.eval(l, i, x);
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
            } else if (t < 2.0) {
              t += 1.0;
              double result = 1.0970809878086267e-02;
              result = -9.5624822442993501e-02 + result * t;
              result = 2.3272870596768197e-01 + result * t;
              result *= t;
              result = -3.8954061556102132e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -1.0565590992210227e-02;
              result = 3.6024896094041699e-02 + result * t;
              result = -3.5470962602601144e-02 + result * t;
              result = -6.5050332141593551e-04 + result * t;
              result = 1.1783132312279444e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 1.5658247612620186e-03;
              result = -6.2374678747992105e-03 + result * t;
              result = 9.2101797262625906e-03 + result * t;
              result = -5.7801042133340313e-03 + result * t;
              result = 1.1209714900938378e-03 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -2.2792209043115266e-06;
              result = 2.5831170248863967e-05 + result * t;
              result = -1.0727533056292918e-04 + result * t;
              result = 1.9115065984159335e-04 + result * t;
              result = -1.2059611051479432e-04 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 1.1048001078675925e-02;
              result = -9.6136991842513314e-02 + result * t;
              result = 2.3353333782781635e-01 + result * t;
              result *= t;
              result = -3.8931211937346905e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -1.0795186408683341e-02;
              result = 3.6439021101597788e-02 + result * t;
              result = -3.5107530506303511e-02 + result * t;
              result = -1.3146362839614374e-03 + result * t;
              result = 1.1677228701768579e-02 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 1.7416121710600238e-03;
              result = -6.7417245331355758e-03 + result * t;
              result = 9.4384143463898066e-03 + result * t;
              result = -5.3933796265084608e-03 + result * t;
              result = 8.9889660441807683e-04 + result * t;
              return innerDeriv * result;
            } else {
              t -= 4.0;
              double result = -5.6181037776129802e-05;
              result = 2.2472415110451921e-04 + result * t;
              result = -3.3708622665677880e-04 + result * t;
              result = 2.2472415110451921e-04 + result * t;
              result = -5.6181037776129802e-05 + result * t;
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
          return weaklyFundamentalNakSplineBasisDeriv1.eval(l, i, x);
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
              double result = 8.2122087516619180e-05;
              result = -1.8267520493011576e-03 + result * t;
              result = 1.5949442453138927e-02 + result * t;
              result = -6.5984965881340382e-02 + result * t;
              result = 1.1418333249373620e-01 + result * t;
              result *= t;
              result = -1.4991360495443518e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -8.4650996320802801e-06;
              result = 1.4417805109770253e-04 + result * t;
              result = -8.7629752889562299e-04 + result * t;
              result = 2.0420575019697803e-03 + result * t;
              result = -2.6227806993463345e-04 + result * t;
              result = -4.4389517542570622e-03 + result * t;
              result = 2.8171385268114225e-03 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 1.2742610611393687e-04;
              result = -2.5803391447930288e-03 + result * t;
              result = 2.0286074427391630e-02 + result * t;
              result = -7.4933690967028954e-02 + result * t;
              result = 1.1527364502034521e-01 + result * t;
              result *= t;
              result = -1.2122116798687088e-01 + result * t;
              return innerDeriv * result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -9.5922242875971716e-05;
              result = 4.7788740194145611e-04 + result * t;
              result = -7.3844300112409775e-04 + result * t;
              result = -1.0534746980828942e-04 + result * t;
              result = 1.4316932555734003e-03 + result * t;
              result = -1.3210622134195524e-03 + result * t;
              result = 3.0603023568060030e-04 + result * t;
              return innerDeriv * result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 1.7348022919897944e-05;
              result = -9.7646055314374209e-05 + result * t;
              result = 2.1216036544360691e-04 + result * t;
              result = -1.9869031240955402e-04 + result * t;
              result = 2.5033215678930413e-05 + result * t;
              result = 8.6413436257438717e-05 + result * t;
              result = -4.5164034032454823e-05 + result * t;
              return innerDeriv * result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -1.0822016402601965e-06;
              result = 6.4420822050134528e-06 + result * t;
              result = -1.5849567329795003e-05 + result * t;
              result = 2.0451054619090325e-05 + result * t;
              result = -1.4315738233363228e-05 + result * t;
              result = 4.9082531085816785e-06 + result * t;
              result = -5.4536145650907532e-07 + result * t;
              return innerDeriv * result;
            } else {
              t -= 6.0;
              double result = 8.5212727579543019e-09;
              result = -5.1127636547725818e-08 + result * t;
              result = 1.2781909136931454e-07 + result * t;
              result = -1.7042545515908605e-07 + result * t;
              result = 1.2781909136931454e-07 + result * t;
              result = -5.1127636547725818e-08 + result * t;
              result = 8.5212727579543019e-09 + result * t;
              return innerDeriv * result;
            }
          }
        }

      default:
        return 0.0;
    }
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "WeaklyFundamentalNakSplineModifiedBasisDeriv1: evalDx not implemented"
              << std::endl;
    return -1;
  }

  inline double getIntegral(LT level, IT index) override {
    std::cerr << "WeaklyFundamentalNakSplineModifiedBasisDeriv1: getIntegral not implemented"
              << std::endl;
    return -1.0;
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override {
    return weaklyFundamentalNakSplineBasisDeriv1.getDegree();
  }

 protected:
  /// Unmodified basis
  WeaklyFundamentalNakSplineBasisDeriv1<LT, IT> weaklyFundamentalNakSplineBasisDeriv1;
};

// default type-def (unsigned int for level and index)
typedef WeaklyFundamentalNakSplineModifiedBasisDeriv1<unsigned int, unsigned int>
    SWeaklyFundamentalNakSplineModifiedBaseDeriv1;

}  // namespace base
}  // namespace sgpp

#endif /* WEAKLY_FUNDAMENTAL_NAK_SPLINE_MODIFIED_BASE_DERIV1_HPP */
