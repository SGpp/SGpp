// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/tools/Distribution.hpp>
#include <sgpp/base/tools/DistributionUniform.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>

namespace sgpp {
namespace base {

/**
 * Not-a-knot polynomial B-spline basis.
 *  "Not-aKn-otter Basis"
 */
template <class LT, class IT>
class NakPBsplineBasis : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakPBsplineBasis() : bsplineBasis(BsplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakPBsplineBasis(size_t degree) : bsplineBasis(BsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 5) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakPBsplineBasis() override {}

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
          // Lagrange level 1
          if (i == 1) {
            return 1.0;
          }
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >=2, 1 < i < 2^l - 1
          return bsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }
          if (i == 1) {
            // l >=2, i = 1
            if ((t < -1.0) || (t > 1.0)) {
              return 0.0;
            } else {
              t += 1.0;
              double result = -1.0 / 2.0;
              result = 1.0 + result * t;
              return result;
            }
          }
        }

      case 3:
        if (l == 1) {
          // Lagrange leve 1
          if (i == 1) {
            // l = 1, i = 1
            // return 4 * x - 4 * x * x;
            return 1.0;
          }
        } else if (l == 2) {
          // Lagrange leve 2
          if (i == 1) {
            // l = 2, i = 1
            return 8 * (x - 0.5) * (x - 0.75);
          } else if (i == 3) {
            // l = 2, i = 3
            return 8 * (x - 0.25) * (x - 0.5);
          }
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 4, 3 < i < 2^l - 3
          return bsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (i == 1) {
            // l >=3, i = 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else {
              t += 1.0;
              double result = -1.0 / 60.0;
              result = 3.0 / 20.0 + result * t;
              result = -9.0 / 20.0 + result * t;
              result = 9.0 / 20.0 + result * t;
              return result;
            }
          } else if (i == 3) {
            // l >= 3, i = 3
            if ((t < -3.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = -1.0 / 20.0;
              result = 3.0 / 20.0 + result * t;
              result = 3.0 / 20.0 + result * t;
              result = 1.0 / 20.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 1.0 / 5.0;
              result = -3.0 / 10.0 + result * t;
              result = -3.0 / 10.0 + result * t;
              result = 1.0 / 2.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0 / 10.0;
              result = 3.0 / 10.0 + result * t;
              result = -3.0 / 10.0 + result * t;
              result = 1.0 / 10.0 + result * t;
              return result;
            }
          }
        }
      case 5:
        if (l == 1) {
          // Lagrange leve 1
          if (i == 1) {
            // l = 1, i = 1
            // return 4 * x - 4 * x * x;
            return 1.0;
          }
        } else if (l == 2) {
          // Lagrange leve 2
          if (i == 1) {
            // l = 2, i = 1
            return 8 * (x - 0.5) * (x - 0.75);
          } else if (i == 3) {
            // l = 2, i = 3
            return 8 * (x - 0.25) * (x - 0.5);
          }
        } else if ((i > 5) && (i < hInv - 5)) {
          // l >= 4, 5 < i < 2^l - 5
          return bsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }
          if (l == 3) {
            if (i == 1) {
              // l=3, i=1
              if ((t < -1.0) || (t > 3.0)) {
                return 0.0;
              } else {
                t += 1.0;
                double result = -1.0 / 6720.0;
                result = 1.0 / 336.0 + result * t;
                result = -1.0 / 42.0 + result * t;
                result = 2.0 / 21.0 + result * t;
                result = -4.0 / 21.0 + result * t;
                result = 16.0 / 105.0 + result * t;
                return result;
              }
            } else if (i == 3) {
              // l=3, i=3
              if ((t < -3.0) || (t > 5.0)) {
                return 0.0;
              } else if (t < 1.0) {
                t += 3.0;
                double result = -149.0 / 277200.0;
                result = 19.0 / 2310.0 + result * t;
                result = -41.0 / 1155.0 + result * t;
                result = -1.0 / 105.0 + result * t;
                result = 499.0 / 2310.0 + result * t;
                result = 2869.0 / 11550.0 + result * t;
                return result;
              } else {
                t -= 1.0;
                double result = 7.0 / 39600.0;
                result = -1.0 / 396.0 + result * t;
                result = 1.0 / 99.0 + result * t;
                result = 1.0 / 99.0 + result * t;
                result = -29.0 / 198.0 + result * t;
                result = 241.0 / 990.0 + result * t;
                return result;
              }
            }
          } else {
            if (i == 1) {
              // l >=4, i = 1
              if ((t < -1.0) || (t > 3.0)) {
                return 0.0;
              } else {
                t += 1.0;
                double result = -1.0 / 6720.0;
                result = 1.0 / 336.0 + result * t;
                result = -1.0 / 42.0 + result * t;
                result = 2.0 / 21.0 + result * t;
                result = -4.0 / 21.0 + result * t;
                result = 16.0 / 105.0 + result * t;
                return result;
              }
            } else if (i == 3) {
              // l >=4, i = 3
              if ((t < -3.0) || (t > 3.0)) {
                return 0.0;
              } else if (t < 1.0) {
                t += 3.0;
                double result = -1.0 / 672.0;
                result = 1.0 / 56.0 + result * t;
                result = -3.0 / 56.0 + result * t;
                result = -3.0 / 56.0 + result * t;
                result = 27.0 / 112.0 + result * t;
                result = 177.0 / 560.0 + result * t;
                return result;
              } else if (t < 2.0) {
                t -= 1.0;
                double result = 13.0 / 3360.0;
                result = -1.0 / 84.0 + result * t;
                result = -1.0 / 168.0 + result * t;
                result = 11.0 / 168.0 + result * t;
                result = -31.0 / 336.0 + result * t;
                result = 71.0 / 1680.0 + result * t;
                return result;
              } else {
                t -= 2.0;
                double result = -1.0 / 672.0;
                result = 5.0 / 672.0 + result * t;
                result = -5.0 / 336.0 + result * t;
                result = 5.0 / 336.0 + result * t;
                result = -5.0 / 672.0 + result * t;
                result = 1.0 / 672.0 + result * t;
                return result;
              }
            } else if (i == 5) {
              // l >=4, i = 5
              if ((t < -5.0) || (t > 3.0)) {
                return 0.0;
              } else if (t < -1.0) {
                t += 5.0;
                double result = -1.0 / 1344.0;
                result = 1.0 / 336.0 + result * t;
                result = 1.0 / 168.0 + result * t;
                result = 1.0 / 168.0 + result * t;
                result = 1.0 / 336.0 + result * t;
                result = 1.0 / 1680.0 + result * t;
                return result;
              } else if (t < 0.0) {
                t += 1.0;
                double result = 121.0 / 6720.0;
                result = -1.0 / 84.0 + result * t;
                result = -11.0 / 168.0 + result * t;
                result = -19.0 / 168.0 + result * t;
                result = 7.0 / 48.0 + result * t;
                result = 821.0 / 1680.0 + result * t;
                return result;
              } else if (t < 1.0) {
                double result = -43.0 / 1344.0;
                result = 5.0 / 64.0 + result * t;
                result = 15.0 / 224.0 + result * t;
                result = -45.0 / 224.0 + result * t;
                result = -15.0 / 64.0 + result * t;
                result = 207.0 / 448.0 + result * t;
                return result;
              } else if (t < 2.0) {
                t -= 1.0;
                double result = 29.0 / 1344.0;
                result = -55.0 / 672.0 + result * t;
                result = 5.0 / 84.0 + result * t;
                result = 25.0 / 168.0 + result * t;
                result = -95.0 / 336.0 + result * t;
                result = 47.0 / 336.0 + result * t;
                return result;
              } else {
                t -= 2.0;
                double result = -1.0 / 192.0;
                result = 5.0 / 192.0 + result * t;
                result = -5.0 / 96.0 + result * t;
                result = 5.0 / 96.0 + result * t;
                result = -5.0 / 192.0 + result * t;
                result = 1.0 / 192.0 + result * t;
                return result;
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
   * @return      value of derivative of B-spline basis function
   */
  inline double evalDx(LT l, IT i, double x) {
    std::cout << "NakPBsplineBasis evalDx not implemented";
    return 0.0;
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @return      integral of basis function
   */
  inline double getIntegral(LT l, IT i) {
    std::cout << "NakPBsplineBasis getIntegral not implemented";
    return -1;
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @return      integral of basis function
   */
  inline double getMean(LT l, IT i, std::shared_ptr<sgpp::base::Distribution> pdf,
                        std::shared_ptr<sgpp::base::DataVector> quadCoordinates,
                        std::shared_ptr<sgpp::base::DataVector> quadWeights) {
    std::cout << "NakPBsplineBasis getMean not implemented";
    return -1;
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const { return bsplineBasis.getDegree(); }

 protected:
  /// B-spline basis for B-spline evaluation
  BsplineBasis<LT, IT> bsplineBasis;

 private:
};  // namespace base
// default type-def (unsigned int for level and index)
typedef NakPBsplineBasis<unsigned int, unsigned int> SNakPBsplineBase;
}  // namespace base
}  // namespace sgpp
