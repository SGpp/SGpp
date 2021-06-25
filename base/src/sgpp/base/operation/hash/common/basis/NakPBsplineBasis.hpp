// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/tools/Distribution.hpp>
#include <sgpp/base/tools/DistributionUniform.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/globaldef.hpp>

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
          return 1.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >=2, 1 < i < 2^l - 1
          return bsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

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
      case 3:
        if (l == 1) {
          // Lagrange leve 1
          // l = 1, i = 1
          return 1.0;
        } else if (l == 2) {
          // Lagrange leve 2
          if (i == 1) {
            // l = 2, i = 1
            return 8 * (x - 0.5) * (x - 0.75);
          } else {
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
          } else {
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
          // l = 1, i = 1
          // return 4 * x - 4 * x * x;
          return 1.0;
        } else if (l == 2) {
          // Lagrange leve 2
          if (i == 1) {
            // l = 2, i = 1
            return 8 * (x - 0.5) * (x - 0.75);
          } else {
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
            } else {
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
            } else {
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
  inline double evalDx(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    double innerDeriv = hInvDbl;
    double t = x * hInvDbl - static_cast<double>(i);
    switch (getDegree()) {
      case 1:
        if (l == 1) {
          // Lagrange level 1
          return 0.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >=2, 1 < i < 2^l - 1
          return bsplineBasis.evalDx(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }
          // l >=2, i = 1
          if ((t < -1.0) || (t > 1.0)) {
            return 0.0;
          } else {
            return -1.0 / 2.0;
          }
        }
      case 3:
        if (l == 1) {
          // Lagrange leve 1
          // l = 1, i = 1
          return 0.0;
        } else if (l == 2) {
          // Lagrange leve 2
          if (i == 1) {
            // l = 2, i = 1
            return 16 * (x - 0.625);
          } else {
            // l = 2, i = 3
            return 16 * (x - 0.375);
          }
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 4, 3 < i < 2^l - 3
          return bsplineBasis.evalDx(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }

          if (i == 1) {
            // l >=3, i = 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else {
              t += 1.0;
              double result = -1.0 / 20.0;
              result = 3.0 / 10.0 + result * t;
              result = -9.0 / 20.0 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 3, i = 3
            if ((t < -3.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = -3.0 / 20.0;
              result = 3.0 / 10.0 + result * t;
              result = 3.0 / 20.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 3.0 / 5.0;
              result = -3.0 / 5.0 + result * t;
              result = -3.0 / 10.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -3.0 / 10.0;
              result = 3.0 / 5.0 + result * t;
              result = -3.0 / 10.0 + result * t;
              return innerDeriv * result;
            }
          }
        }
      case 5:
        if (l == 1) {
          // Lagrange leve 1
          // l = 1, i = 1
          return 0.0;
        } else if (l == 2) {
          // Lagrange leve 2
          if (i == 1) {
            // l = 2, i = 1
            return 16 * (x - 0.625);
          } else {
            // l = 2, i = 3
            return 16 * (x - 0.375);
          }
        } else if ((i > 5) && (i < hInv - 5)) {
          // l >= 4, 5 < i < 2^l - 5
          return bsplineBasis.evalDx(l, i, x);
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
                double result = -1.0 / 1344.0;
                result = 1.0 / 84.0 + result * t;
                result = -1.0 / 14.0 + result * t;
                result = 4.0 / 21.0 + result * t;
                result = -4.0 / 21.0 + result * t;
                return result;
              }
            } else {
              // l=3, i=3
              if ((t < -3.0) || (t > 5.0)) {
                return 0.0;
              } else if (t < 1.0) {
                t += 3.0;
                double result = -149.0 / 55440.0;
                result = 38.0 / 1155.0 + result * t;
                result = -41.0 / 385.0 + result * t;
                result = -2.0 / 105.0 + result * t;
                result = 499.0 / 2310.0 + result * t;
                return result;
              } else {
                t -= 1.0;
                double result = 7.0 / 7920.0;
                result = -1.0 / 99.0 + result * t;
                result = 1.0 / 33.0 + result * t;
                result = 2.0 / 99.0 + result * t;
                result = -29.0 / 198.0 + result * t;
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
                double result = -1.0 / 1344.0;
                result = 1.0 / 84.0 + result * t;
                result = -1.0 / 14.0 + result * t;
                result = 4.0 / 21.0 + result * t;
                result = -4.0 / 21.0 + result * t;
                return result;
              }
            } else if (i == 3) {
              // l >=4, i = 3
              if ((t < -3.0) || (t > 3.0)) {
                return 0.0;
              } else if (t < 1.0) {
                t += 3.0;
                double result = -5.0 / 672.0;
                result = 1.0 / 14.0 + result * t;
                result = -9.0 / 56.0 + result * t;
                result = -3.0 / 28.0 + result * t;
                result = 27.0 / 112.0 + result * t;
                return result;
              } else if (t < 2.0) {
                t -= 1.0;
                double result = 13.0 / 672.0;
                result = -1.0 / 21.0 + result * t;
                result = -1.0 / 56.0 + result * t;
                result = 11.0 / 84.0 + result * t;
                result = -31.0 / 336.0 + result * t;
                return result;
              } else {
                t -= 2.0;
                double result = -5.0 / 672.0;
                result = 5.0 / 168.0 + result * t;
                result = -5.0 / 112.0 + result * t;
                result = 5.0 / 168.0 + result * t;
                result = -5.0 / 672.0 + result * t;
                return result;
              }

            } else {
              // l >=4, i = 5
              if ((t < -5.0) || (t > 3.0)) {
                return 0.0;
              } else if (t < -1.0) {
                t += 5.0;
                double result = -5.0 / 1344.0;
                result = 1.0 / 84.0 + result * t;
                result = 1.0 / 56.0 + result * t;
                result = 1.0 / 84.0 + result * t;
                result = 1.0 / 336.0 + result * t;
                return result;
              } else if (t < 0.0) {
                t += 1.0;
                double result = 121.0 / 1344.0;
                result = -1.0 / 21.0 + result * t;
                result = -11.0 / 56.0 + result * t;
                result = -19.0 / 84.0 + result * t;
                result = 7.0 / 48.0 + result * t;
                return result;
              } else if (t < 1.0) {
                double result = -215.0 / 1344.0;
                result = 5.0 / 16.0 + result * t;
                result = 45.0 / 224.0 + result * t;
                result = -45.0 / 112.0 + result * t;
                result = -15.0 / 64.0 + result * t;
                return result;
              } else if (t < 2.0) {
                t -= 1.0;
                double result = 145.0 / 1344.0;
                result = -55.0 / 168.0 + result * t;
                result = 5.0 / 28.0 + result * t;
                result = 25.0 / 84.0 + result * t;
                result = -95.0 / 336.0 + result * t;
                return result;
              } else {
                t -= 2.0;
                double result = -5.0 / 192.0;
                result = 5.0 / 48.0 + result * t;
                result = -5.0 / 32.0 + result * t;
                result = 5.0 / 48.0 + result * t;
                result = -5.0 / 192.0 + result * t;
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
   * @return      integral of basis function
   */
  inline double getIntegral(LT l, IT i) override {
    size_t quadOrder = getDegree() + 1;
    auto pdf_uniform = std::make_shared<sgpp::base::DistributionUniform>(0, 1);
    base::DataVector temp_quadCoordinates, temp_quadWeights;
    base::GaussLegendreQuadRule1D gauss;
    gauss.getLevelPointsAndWeightsNormalized(quadOrder, temp_quadCoordinates, temp_quadWeights);
    auto quadCoordinates = std::make_shared<sgpp::base::DataVector>(temp_quadCoordinates);
    auto quadWeights = std::make_shared<sgpp::base::DataVector>(temp_quadWeights);
    return getMean(l, i, pdf_uniform, quadCoordinates, quadWeights);
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param pdf   probability density function
   * @param quadCoordinates coordinates of the quadrature rule to be used
   * @param quadWeights weights of the quadrature rule to be used
   * @return      integral of basis function
   */
  inline double getMean(LT l, IT i, std::shared_ptr<sgpp::base::Distribution> pdf,
                        std::shared_ptr<sgpp::base::DataVector> quadCoordinates,
                        std::shared_ptr<sgpp::base::DataVector> quadWeights) {
    size_t degree = getDegree();
    if ((degree != 1) && (degree != 3) && (degree != 5)) {
      throw std::runtime_error("NakPBsplineBasis: only B spline degrees 1, 3 and 5 are supported.");
    }

    const size_t pp1h = (degree + 1) >> 1;  //  =|_(p+1)/2_|
    const size_t hInv = 1 << l;             // = 2^lid
    const double hik = 1.0 / static_cast<double>(hInv);
    double offset = (i - static_cast<double>(pp1h)) * hik;

    if (degree == 3) {
      if (i == 3) offset -= hik;
    } else if (degree == 5) {
      if (i == 5) offset -= 2 * hik;
    }
    // start and stop identify the segments on which the spline is nonzero
    size_t start = 0, stop = 0;
    start = ((i > pp1h) ? 0 : (pp1h - i));
    stop = std::min(degree, hInv + pp1h - i - 1);
    // nak special cases
    // TODO (rehmemk) Copied from nakBsplineExtended. Verify that this holds here too
    if (degree == 3) {
      if ((i == 3) || (i == hInv - 3)) stop += 1;
    } else if (degree == 5) {
      if ((i == 5) || (i == hInv - 5)) stop += 2;
    }
    if (l == 2) {
      start = 1;
      stop = 4;
      offset = -0.25;
    }
    if ((degree == 5) && (l == 3)) {
      start = 1;
      stop = 8;
      offset = -0.125;
    }

    double temp_res = basisMean(l, i, start, stop, offset, hik, quadCoordinates, quadWeights, pdf);
    double integral = temp_res * hik;

    return integral;
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override { return bsplineBasis.getDegree(); }

 protected:
  /// B-spline basis for B-spline evaluation
  BsplineBasis<LT, IT> bsplineBasis;

 private:
  double basisMean(LT l, IT i, size_t start, size_t stop, double offset, double hik,
                   std::shared_ptr<sgpp::base::DataVector> quadCoordinates,
                   std::shared_ptr<sgpp::base::DataVector> quadWeights,
                   std::shared_ptr<sgpp::base::Distribution> pdf) {
    sgpp::base::DataVector bounds = pdf->getBounds();
    double left = bounds[0];
    double right = bounds[1];

    double temp_res = 0.0;
    // loop over the segments the B-spline is defined on
    for (size_t n = start; n <= stop; n++) {
      // loop over quadrature points
      for (size_t c = 0; c < quadCoordinates->getSize(); c++) {
        // transform  the quadrature points to the segment on which the Bspline is
        // evaluated and the support of the pdf
        double x = offset + hik * (quadCoordinates->get(c) + static_cast<double>(n));
        double scaledX = left + (right - left) * x;
        temp_res += quadWeights->get(c) * this->eval(l, i, x) * pdf->eval(scaledX);
      }
    }
    return temp_res * (right - left);
  }
};
// default type-def (unsigned int for level and index)
typedef NakPBsplineBasis<unsigned int, unsigned int> SNakPBsplineBase;
}  // namespace base
}  // namespace sgpp
