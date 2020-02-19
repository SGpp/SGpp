// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WEAKLY_FUNDAMENTAL_NAK_SPLINE_MODIFIED_BASE_HPP
#define WEAKLY_FUNDAMENTAL_NAK_SPLINE_MODIFIED_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * Modified weakly fundamental not-a-knot spline basis.
 */
template <class LT, class IT>
class WeaklyFundamentalNakSplineModifiedBasis : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  WeaklyFundamentalNakSplineModifiedBasis()
      : weaklyFundamentalNakSplineBasis(WeaklyFundamentalNakSplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    Spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit WeaklyFundamentalNakSplineModifiedBasis(size_t degree)
      : weaklyFundamentalNakSplineBasis(WeaklyFundamentalNakSplineBasis<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("Unsupported spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~WeaklyFundamentalNakSplineModifiedBasis() override {}

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
          return weaklyFundamentalNakSplineBasis.eval(l, i, x);
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
              double result = 4.3269230769230768e-02;
              result *= t;
              result = -8.0769230769230771e-01 + result * t;
              result = 1.2692307692307692e+00 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -5.7692307692307696e-02;
              result = 2.5961538461538464e-01 + result * t;
              result = -2.8846153846153844e-01 + result * t;
              result *= t;
              return result;
            }
          } else {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 5.3571428571428568e-02;
              result *= t;
              result = -8.5714285714285710e-01 + result * t;
              result = 1.2857142857142858e+00 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = -1.2500000000000000e-01;
              result = 3.2142857142857145e-01 + result * t;
              result = -2.1428571428571427e-01 + result * t;
              result *= t;
              return result;
            } else {
              t -= 2.0;
              double result = 1.7857142857142856e-02;
              result = -5.3571428571428568e-02 + result * t;
              result = 5.3571428571428568e-02 + result * t;
              result = -1.7857142857142856e-02 + result * t;
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
          return weaklyFundamentalNakSplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -3.0 / 70.0;
            result = 19.0 / 210.0 + result * t;
            result = 37.0 / 70.0 + result * t;
            result = -331.0 / 210.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else if (l == 3) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 2.1941619756172532e-03;
              result = -2.3906205610748375e-02 + result * t;
              result = 7.7576235322560652e-02 + result * t;
              result *= t;
              result = -3.8954061556102132e-01 + result * t;
              result = 4.7075745509377931e-01 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -2.1131181984420455e-03;
              result = 9.0062240235104247e-03 + result * t;
              result = -1.1823654200867048e-02 + result * t;
              result = -3.2525166070796775e-04 + result * t;
              result = 1.1783132312279444e-02 + result * t;
              result = -6.5273322757728076e-03 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 3.1316495225240372e-04;
              result = -1.5593669686998026e-03 + result * t;
              result = 3.0700599087541969e-03 + result * t;
              result = -2.8900521066670157e-03 + result * t;
              result = 1.1209714900938378e-03 + result * t;
              result *= t;
              return result;
            } else {
              t -= 4.0;
              double result = -4.5584418086230533e-07;
              result = 6.4577925622159919e-06 + result * t;
              result = -3.5758443520976394e-05 + result * t;
              result = 9.5575329920796675e-05 + result * t;
              result = -1.2059611051479432e-04 + result * t;
              result = 5.4777275733620351e-05 + result * t;
              return result;
            }
          } else {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 2.2096002157351849e-03;
              result = -2.4034247960628329e-02 + result * t;
              result = 7.7844445942605450e-02 + result * t;
              result *= t;
              result = -3.8931211937346905e-01 + result * t;
              result = 4.6970943167262186e-01 + result * t;
              return result;
            } else if (t < 3.0) {
              t -= 2.0;
              double result = -2.1590372817366682e-03;
              result = 9.1097552753994469e-03 + result * t;
              result = -1.1702510168767838e-02 + result * t;
              result = -6.5731814198071868e-04 + result * t;
              result = 1.1677228701768579e-02 + result * t;
              result = -6.2681183846828017e-03 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = 3.4832243421200476e-04;
              result = -1.6854311332838939e-03 + result * t;
              result = 3.1461381154632689e-03 + result * t;
              result = -2.6966898132542304e-03 + result * t;
              result = 8.9889660441807683e-04 + result * t;
              result *= t;
              return result;
            } else {
              t -= 4.0;
              double result = -1.1236207555225960e-05;
              result = 5.6181037776129802e-05 + result * t;
              result = -1.1236207555225960e-04 + result * t;
              result = 1.1236207555225960e-04 + result * t;
              result = -5.6181037776129802e-05 + result * t;
              result = 1.1236207555225960e-05 + result * t;
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
          return weaklyFundamentalNakSplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1
            double result = -3.0 / 70.0;
            result = 19.0 / 210.0 + result * t;
            result = 37.0 / 70.0 + result * t;
            result = -331.0 / 210.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else if (l == 3) {
            // l = 3, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 1.1731726788088454e-05;
              result = -3.0445867488352625e-04 + result * t;
              result = 3.1898884906277855e-03 + result * t;
              result = -1.6496241470335096e-02 + result * t;
              result = 3.8061110831245401e-02 + result * t;
              result *= t;
              result = -1.4991360495443518e-01 + result * t;
              result = 1.7518544924784990e-01 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.2092999474400401e-06;
              result = 2.4029675182950422e-05 + result * t;
              result = -1.7525950577912462e-04 + result * t;
              result = 5.1051437549244508e-04 + result * t;
              result = -8.7426023311544480e-05 + result * t;
              result = -2.2194758771285311e-03 + result * t;
              result = 2.8171385268114225e-03 + result * t;
              result *= t;
              return result;
            }
          } else {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 7.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 1.8203729444848127e-05;
              result = -4.3005652413217148e-04 + result * t;
              result = 4.0572148854783260e-03 + result * t;
              result = -1.8733422741757239e-02 + result * t;
              result = 3.8424548340115064e-02 + result * t;
              result *= t;
              result = -1.2122116798687088e-01 + result * t;
              result = 1.3014337696114900e-01 + result * t;
              return result;
            } else if (t < 4.0) {
              t -= 3.0;
              double result = -1.3703177553710245e-05;
              result = 7.9647900323576019e-05 + result * t;
              result = -1.4768860022481957e-04 + result * t;
              result = -2.6336867452072355e-05 + result * t;
              result = 4.7723108519113344e-04 + result * t;
              result = -6.6053110670977622e-04 + result * t;
              result = 3.0603023568060030e-04 + result * t;
              result *= t;
              return result;
            } else if (t < 5.0) {
              t -= 4.0;
              double result = 2.4782889885568491e-06;
              result = -1.6274342552395704e-05 + result * t;
              result = 4.2432073088721378e-05 + result * t;
              result = -4.9672578102388505e-05 + result * t;
              result = 8.3444052263101371e-06 + result * t;
              result = 4.3206718128719359e-05 + result * t;
              result = -4.5164034032454823e-05 + result * t;
              result = 1.4649469254931299e-05 + result * t;
              return result;
            } else if (t < 6.0) {
              t -= 5.0;
              double result = -1.5460023432288520e-07;
              result = 1.0736803675022421e-06 + result * t;
              result = -3.1699134659590005e-06 + result * t;
              result = 5.1127636547725814e-06 + result * t;
              result = -4.7719127444544091e-06 + result * t;
              result = 2.4541265542908393e-06 + result * t;
              result = -5.4536145650907532e-07 + result * t;
              result *= t;
              return result;
            } else {
              t -= 6.0;
              double result = 1.2173246797077576e-09;
              result = -8.5212727579543019e-09 + result * t;
              result = 2.5563818273862909e-08 + result * t;
              result = -4.2606363789771513e-08 + result * t;
              result = 4.2606363789771513e-08 + result * t;
              result = -2.5563818273862909e-08 + result * t;
              result = 8.5212727579543019e-09 + result * t;
              result = -1.2173246797077576e-09 + result * t;
              return result;
            }
          }
        }

      default:
        return 0.0;
    }
  }

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "WeaklyFundamentalNakSplineModifiedBasis: evalDx not implemented" << std::endl;
    return -1;
  }

  inline double getIntegral(LT level, IT index) override {
    std::cerr << "WeaklyFundamentalNakSplineModifiedBasis: getIntegral not implemented"
              << std::endl;
    return -1.0;
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override { return weaklyFundamentalNakSplineBasis.getDegree(); }

 protected:
  /// Unmodified basis
  WeaklyFundamentalNakSplineBasis<LT, IT> weaklyFundamentalNakSplineBasis;
};

// default type-def (unsigned int for level and index)
typedef WeaklyFundamentalNakSplineModifiedBasis<unsigned int, unsigned int>
    SWeaklyFundamentalNakSplineModifiedBase;

}  // namespace base
}  // namespace sgpp

#endif /* WEAKLY_FUNDAMENTAL_NAK_SPLINE_MODIFIED_BASE_HPP */
