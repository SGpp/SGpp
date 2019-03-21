// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/Distribution.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
// #include <sgpp/combigrid/GeneralFunction.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include "NakBsplineBasis.hpp"

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
  NakBsplineModifiedBasis() : notAKnotBsplineBasis(NakBsplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineModifiedBasis(size_t degree)
      : notAKnotBsplineBasis(NakBsplineBasis<LT, IT>(degree)) {
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

          // l >= 3, i = 1
          return std::max(1.0 - t, 0.0);
        }

      case 3:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return notAKnotBsplineBasis.eval(l, i, x);
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
          return notAKnotBsplineBasis.eval(l, i, x);
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
                double result = 37.0 / 22680.0;
                result = -1.0 / 54.0 + result * t;
                result = 4.0 / 63.0 + result * t;
                result *= t;
                result = -8.0 / 21.0 + result * t;
                result = 8.0 / 15.0 + result * t;
                return result;
              } else {
                t -= 2.0;
                double result = -1.0 / 840.0;
                result = 1.0 / 168.0 + result * t;
                result = -1.0 / 84.0 + result * t;
                result = 1.0 / 84.0 + result * t;
                result = -1.0 / 168.0 + result * t;
                result = 1.0 / 840.0 + result * t;
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

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double evalDx(LT l, IT i, double x) {
    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    double innerDeriv = hInvDbl;
    double t = x * hInvDbl - static_cast<double>(i);

    switch (getDegree()) {
      case 1:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return notAKnotBsplineBasis.evalDx(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
            innerDeriv *= -1.0;
          }

          // l >= 3, i = 1
          if (t > 1.0) {
            return 0.0;
          } else {
            return -innerDeriv;
          }
        }

      case 3:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if ((i > 1) && (i < hInv - 1)) {
          // l >= 3, 1 < i < 2^l - 1
          return notAKnotBsplineBasis.evalDx(l, i, x);
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
              double result = 7.4999999999999997e-02;
              result *= t;
              result = -5.9999999999999998e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -7.4999999999999997e-02;
              result = 2.9999999999999999e-01 + result * t;
              result = -2.9999999999999999e-01 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 1.2500000000000000e-01;
              result *= t;
              result = -7.5000000000000000e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 1.0;
              double result = -2.5000000000000000e-01;
              result = 5.0000000000000000e-01 + result * t;
              result = -2.5000000000000000e-01 + result * t;
              return innerDeriv * result;
            }
          }
        }

      case 5:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 3, 3 < i < 2^l - 3
          return notAKnotBsplineBasis.evalDx(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (l == 2) {
            // l = 2, i = 1 (Lagrange)
            double result = -2.0 / 23.0;
            result = 6.0 / 23.0 + result * t;
            result = 18.0 / 23.0 + result * t;
            result = -67.0 / 46.0 + result * t;
            return result;
          } else if ((l == 3) && (i == 3)) {
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 1.3778659611992945e-03;
              result = 2.7513227513227514e-02 + result * t;
              result = -1.6825396825396827e-01 + result * t;
              result *= t;
              result = 3.4973544973544973e-01 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -3.6276455026455025e-02;
              result = 4.4047619047619051e-02 + result * t;
              result = 1.5376984126984128e-01 + result * t;
              result = -1.1785714285714285e-01 + result * t;
              result = -3.1008597883597883e-01 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 2.8736772486772488e-02;
              result = -1.0105820105820106e-01 + result * t;
              result = 6.8253968253968247e-02 + result * t;
              result = 1.7671957671957672e-01 + result * t;
              result = -2.6640211640211642e-01 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.1574074074074073e-03;
              result = 1.3888888888888888e-02 + result * t;
              result = -6.2500000000000000e-02 + result * t;
              result = 1.2500000000000000e-01 + result * t;
              result = -9.3750000000000000e-02 + result * t;
              return result;
            }
          } else {
            if (i == 1) {
              // l >= 3, i = 1
              if ((t < -1.0) || (t > 3.0)) {
                return 0.0;
              } else if (t < 2.0) {
                t += 1.0;
                double result = 8.1569664902998232e-03;
                result = -7.4074074074074070e-02 + result * t;
                result = 1.9047619047619047e-01 + result * t;
                result *= t;
                result = -3.8095238095238093e-01 + result * t;
                return result;
              } else {
                t -= 2.0;
                double result = -5.9523809523809521e-03;
                result = 2.3809523809523808e-02 + result * t;
                result = -3.5714285714285712e-02 + result * t;
                result = 2.3809523809523808e-02 + result * t;
                result = -5.9523809523809521e-03 + result * t;
                return result;
              }
            } else {
              // l >= 4, i = 3
              if ((t < -3.0) || (t > 3.0)) {
                return 0.0;
              } else if (t < 0.0) {
                t += 3.0;
                double result = 3.9682539682539680e-03;
                result = 1.9047619047619049e-02 + result * t;
                result = -1.7142857142857143e-01 + result * t;
                result *= t;
                result = 3.4285714285714286e-01 + result * t;
                return result;
              } else if (t < 1.0) {
                double result = -6.2698412698412698e-02;
                result = 6.6666666666666666e-02 + result * t;
                result = 2.1428571428571427e-01 + result * t;
                result = -8.5714285714285715e-02 + result * t;
                result = -3.6428571428571427e-01 + result * t;
                return result;
              } else if (t < 2.0) {
                t -= 1.0;
                double result = 6.5873015873015875e-02;
                result = -1.8412698412698414e-01 + result * t;
                result = 3.8095238095238099e-02 + result * t;
                result = 2.9206349206349208e-01 + result * t;
                result = -2.3174603174603176e-01 + result * t;
                return result;
              } else {
                t -= 2.0;
                double result = -1.9841269841269840e-02;
                result = 7.9365079365079361e-02 + result * t;
                result = -1.1904761904761904e-01 + result * t;
                result = 7.9365079365079361e-02 + result * t;
                result = -1.9841269841269840e-02 + result * t;
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
  inline double getIntegral(LT l, IT i) {
    size_t degree = getDegree();
    if ((degree != 1) && (degree != 3) && (degree != 5)) {
      throw std::runtime_error(
          "NakBsplineModified: only B spline degrees 1, 3 and 5 are "
          "supported.");
    }
    size_t quadOrder = degree + 1;
    base::DataVector quadCoordinates, quadWeights;
    base::GaussLegendreQuadRule1D gauss;

    const size_t pp1h = (degree + 1) >> 1;  //  =|_(p+1)/2_|
    const size_t hInv = 1 << l;             // = 2^lid
    const double hik = 1.0 / static_cast<double>(hInv);
    double offset = (i - static_cast<double>(pp1h)) * hik;

    if (degree == 3) {
      if (i == 3) offset -= hik;
    } else if (degree == 5) {
      if (i == 5) offset -= 2 * hik;
    }
    gauss.getLevelPointsAndWeightsNormalized(quadOrder, quadCoordinates, quadWeights);
    // start and stop identify the segments on which the spline is nonzero
    size_t start = 0, stop = 0;
    double scaling = hik;
    start = ((i > pp1h) ? 0 : (pp1h - i));
    stop = std::min(degree, hInv + pp1h - i - 1);
    // nak special cases
    if (degree == 3) {
      if ((i == 3) || (i == hInv - 3)) stop += 1;
    } else if (degree == 5) {
      if ((i == 5) || (i == hInv - 5)) stop += 2;
    }
    if (l == 2) {
      start = 1;
      stop = 4;
      offset = -0.25;
      scaling = 0.25;
    }
    if ((degree == 5) && (l == 3)) {
      start = 1;
      stop = 8;
      offset = -0.125;
      scaling = 0.125;
    }

    double temp_res =
        integrateBspline(l, i, start, stop, offset, scaling, quadCoordinates, quadWeights);
    double integral = temp_res * scaling;

    return integral;
  }

  /**
   * Calculates the mean \int b_i(x) \rho(x) dx of a basis function b_i w.r.t. the probability
   * density function \rho
   *
   * @param l     		level of basis function
   * @param i     		index of basis function
   * @param pdf   		probability density function
   * @param quadOrder	order of the Gauss Legendre quadrature
   * @return      		mean of basis function
   */
  inline double getMean(LT l, IT i, std::shared_ptr<sgpp::base::Distribution> pdf,
                        size_t quadOrder) {
    size_t degree = getDegree();
    if ((degree != 1) && (degree != 3) && (degree != 5)) {
      throw std::runtime_error(
          "NakBsplineModified: only B spline degrees 1, 3 and 5 are "
          "supported.");
    }

    base::DataVector quadCoordinates, quadWeights;
    base::GaussLegendreQuadRule1D gauss;

    const size_t pp1h = (degree + 1) >> 1;  //  =|_(p+1)/2_|
    const size_t hInv = 1 << l;             // = 2^lid
    const double hik = 1.0 / static_cast<double>(hInv);
    double offset = (i - static_cast<double>(pp1h)) * hik;

    if (degree == 3) {
      if (i == 3) offset -= hik;
    } else if (degree == 5) {
      if (i == 5) offset -= 2 * hik;
    }
    gauss.getLevelPointsAndWeightsNormalized(quadOrder, quadCoordinates, quadWeights);
    // start and stop identify the segments on which the spline is nonzero
    size_t start = 0, stop = 0;
    start = ((i > pp1h) ? 0 : (pp1h - i));
    stop = std::min(degree, hInv + pp1h - i - 1);
    // nak special cases
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
    double mean = temp_res * hik;

    return mean;
  }

  double basisMean(LT l, IT i, size_t start, size_t stop, double offset, double hik,
                   base::DataVector quadCoordinates, base::DataVector quadWeights,
                   std::shared_ptr<sgpp::base::Distribution> pdf) {
    double temp_res = 0.0;
    // loop over the segments the B-spline is defined on
    for (size_t n = start; n <= stop; n++) {
      // loop over quadrature points
      for (size_t c = 0; c < quadCoordinates.getSize(); c++) {
        // transform  the quadrature points to the segment on which the Bspline is
        // evaluated
        const double x = offset + hik * (quadCoordinates[c] + static_cast<double>(n));
        temp_res += quadWeights[c] * this->eval(l, i, x) * pdf->eval(x);
      }
    }
    return temp_res;
  }

  //#ifdef SG_COMBIGRID
  /**
   * @param l     				level of basis function
   * @param i     				index of basis function
   * @param weightfunction		weightfunction (usually a probability density function)
   * @param lbound				left boundary of the definition range
   * @param rbound 				right boundary of the definition range
   * @param numAdditionalPoints	number of additional points for the integration of the weighted
   * 							spline
   * @param incrementQuadraturePoints	increment to numAdditionalPoints in each loop cycle
   * @return      				integral of basis function times weight function
   */
  //  inline double getWeightedIntegral(LT l, IT i, sgpp::combigrid::SingleFunction weightfunction,
  //                                    double lbound, double rbound, size_t numAdditionalPoints =
  //                                    0, size_t incrementQuadraturePoints = 10) {
  //    size_t degree = getDegree();
  //    if ((degree != 1) && (degree != 3) && (degree != 5)) {
  //      throw std::runtime_error(
  //          "OperationMatrixLTwoDotNakBsplineBoundary: only B spline degrees 1, 3 and 5 are "
  //          "supported.");
  //    }
  //    size_t quadOrder = degree + 1 + numAdditionalPoints;
  //    base::DataVector quadCoordinates, quadWeights;
  //    base::GaussLegendreQuadRule1D gauss;
  //
  //    const size_t pp1h = (degree + 1) >> 1;  //  =|_(p+1)/2_|
  //    const size_t hInv = 1 << l;             // = 2^lid
  //    const double hik = 1.0 / static_cast<double>(hInv);
  //    double offset = (i - static_cast<double>(pp1h)) * hik;
  //
  //    if (degree == 3) {
  //      if (i == 3) offset -= hik;
  //    } else if (degree == 5) {
  //      if (i == 5) offset -= 2 * hik;
  //    }
  //    gauss.getLevelPointsAndWeightsNormalized(quadOrder, quadCoordinates, quadWeights);
  //    // start and stop identify the segments on which the spline is nonzero
  //    size_t start = 0, stop = 0;
  //    double scaling = hik;
  //    start = ((i > pp1h) ? 0 : (pp1h - i));
  //    stop = std::min(degree, hInv + pp1h - i - 1);
  //    // nak special cases
  //    if (degree == 3) {
  //      if ((i == 3) || (i == hInv - 3)) stop += 1;
  //    } else if (degree == 5) {
  //      if ((i == 5) || (i == hInv - 5)) stop += 2;
  //    }
  //    if (l == 2) {
  //      start = 1;
  //      stop = 4;
  //      offset = -0.25;
  //      scaling = 0.25;
  //    }
  //    if ((degree == 5) && (l == 3)) {
  //      start = 1;
  //      stop = 8;
  //      offset = -0.125;
  //      scaling = 0.125;
  //    }
  //
  //    double temp_res = integrateWeightedBspline(l, i, start, stop, offset, scaling,
  //    quadCoordinates,
  //                                               quadWeights, weightfunction);
  //    double width = rbound - lbound;
  //    temp_res *= width;
  //
  //    double tol = 1e-14;
  //    double err = 1e14;
  //    while (err > tol) {
  //      numAdditionalPoints += incrementQuadraturePoints;
  //      quadOrder = degree + 1 + numAdditionalPoints;
  //      if (quadOrder > 480) {
  //        break;
  //      }
  //      gauss.getLevelPointsAndWeightsNormalized(quadOrder, quadCoordinates, quadWeights);
  //      double finer_temp_res = integrateWeightedBspline(
  //          l, i, start, stop, offset, scaling, quadCoordinates, quadWeights, weightfunction);
  //      finer_temp_res *= width;
  //      err = fabs(temp_res - finer_temp_res);
  //      temp_res = finer_temp_res;
  //    }
  //
  //    double integral = temp_res * scaling;
  //    return integral;
  //  }
  //#endif

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const { return notAKnotBsplineBasis.getDegree(); }

 protected:
  /// B-spline basis for B-spline evaluation
  NakBsplineBasis<LT, IT> notAKnotBsplineBasis;

 private:
  /**
   * integrate one basis function
   *
   * @param l				level
   * @param i				index
   * @param start			index for the supports first segment (usually 0)
   * @param stop			index for the supports last segment (usually degree)
   * @param offset			left point of the support
   * @param scaling			size of one support segment
   * @param quadCoordinates	the quadrature points
   * @param quadWeights		the quadrature weights
   *
   * @return integral
   */
  double integrateBspline(LT l, IT i, size_t start, size_t stop, double offset, double scaling,
                          base::DataVector quadCoordinates, base::DataVector quadWeights) {
    double temp_res = 0.0;
    // loop over the segments the B-spline is defined on
    for (size_t n = start; n <= stop; n++) {
      // loop over quadrature points
      for (size_t c = 0; c < quadCoordinates.getSize(); c++) {
        // transform  the quadrature points to the segment on which the Bspline is
        // evaluated
        const double x = offset + scaling * (quadCoordinates[c] + static_cast<double>(n));
        temp_res += quadWeights[c] * this->eval(l, i, x);
      }
    }
    return temp_res;
  }

  //#ifdef SG_COMBIGRID
  //  double integrateWeightedBspline(LT l, IT i, size_t start, size_t stop, double offset,
  //                                  double scaling, base::DataVector quadCoordinates,
  //                                  base::DataVector quadWeights,
  //                                  sgpp::combigrid::SingleFunction weightfunction) {
  //    double temp_res = 0.0;
  //    for (size_t n = start; n <= stop; n++) {
  //      for (size_t c = 0; c < quadCoordinates.getSize(); c++) {
  //        // transform  the quadrature points to the segment on which the Bspline is
  //        // evaluated
  //        const double x = offset + scaling * (quadCoordinates[c] + static_cast<double>(n));
  //        if (this->eval(l, i, x) == 0) {
  //        }
  //        temp_res += quadWeights[c] * this->eval(l, i, x) * weightfunction(x);
  //      }
  //    }
  //    return temp_res;
  //  }
  //#endif
};  // namespace base

// default type-def (unsigned int for level and index)
typedef NakBsplineModifiedBasis<unsigned int, unsigned int> SNakBsplineModifiedBase;

}  // namespace base
}  // namespace sgpp
