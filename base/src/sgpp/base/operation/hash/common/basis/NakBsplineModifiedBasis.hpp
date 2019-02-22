// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
//#include <sgpp/combigrid/GeneralFunction.hpp>

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
            double result = -3.0 / 70.0;
            result = 19.0 / 210.0 + result * t;
            result = 37.0 / 70.0 + result * t;
            result = -331.0 / 210.0 + result * t;
            result = 1.0 + result * t;
            return result;
          } else {
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
          }
        }

      case 7:
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
              double result = 31.0 / 4942080.0;
              result = -7.0 / 41184.0 + result * t;
              result = 49.0 / 25740.0 + result * t;
              result = -14.0 / 1287.0 + result * t;
              result = 112.0 / 3861.0 + result * t;
              result *= t;
              result = -3584.0 / 19305.0 + result * t;
              result = 2048.0 / 6435.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0 / 4942080.0;
              result = 7.0 / 1235520.0 + result * t;
              result = -7.0 / 102960.0 + result * t;
              result = 7.0 / 15444.0 + result * t;
              result = -7.0 / 3861.0 + result * t;
              result = 28.0 / 6435.0 + result * t;
              result = -112.0 / 19305.0 + result * t;
              result = 64.0 / 19305.0 + result * t;
              return result;
            }
          } else {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 2101.0 / 154828800.0;
              result = -41.0 / 122880.0 + result * t;
              result = 61.0 / 18432.0 + result * t;
              result = -25.0 / 1536.0 + result * t;
              result = 125.0 / 3456.0 + result * t;
              result *= t;
              result = -125.0 / 864.0 + result * t;
              result = 125.0 / 672.0 + result * t;
              return result;
            } else {
              t -= 3.0;
              double result = -1.0 / 151200.0;
              result = 1.0 / 21600.0 + result * t;
              result = -1.0 / 7200.0 + result * t;
              result = 1.0 / 4320.0 + result * t;
              result = -1.0 / 4320.0 + result * t;
              result = 1.0 / 7200.0 + result * t;
              result = -1.0 / 21600.0 + result * t;
              result = 1.0 / 151200.0 + result * t;
              return result;
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
            double result = -6.0 / 35.0;
            result = 19.0 / 70.0 + result * t;
            result = 37.0 / 35.0 + result * t;
            result = -331.0 / 210.0 + result * t;
            return innerDeriv * result;
          } else {
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
              return innerDeriv * result;
            } else {
              t -= 2.0;
              double result = -5.9523809523809521e-03;
              result = 2.3809523809523808e-02 + result * t;
              result = -3.5714285714285712e-02 + result * t;
              result = 2.3809523809523808e-02 + result * t;
              result = -5.9523809523809521e-03 + result * t;
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
          return notAKnotBsplineBasis.evalDx(l, i, x);
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
              double result = 4.3908637658637657e-05;
              result = -1.0198135198135198e-03 + result * t;
              result = 9.5182595182595180e-03 + result * t;
              result = -4.3512043512043512e-02 + result * t;
              result = 8.7024087024087024e-02 + result * t;
              result *= t;
              result = -1.8565138565138564e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -1.4164076664076663e-06;
              result = 3.3993783993783994e-05 + result * t;
              result = -3.3993783993783994e-04 + result * t;
              result = 1.8130018130018131e-03 + result * t;
              result = -5.4390054390054390e-03 + result * t;
              result = 8.7024087024087024e-03 + result * t;
              result = -5.8016058016058013e-03 + result * t;
              return innerDeriv * result;
            }
          } else {
            // l >= 4, i = 1
            if ((t < -1.0) || (t > 4.0)) {
              return 0.0;
            } else if (t < 3.0) {
              t += 1.0;
              double result = 9.4988787615740740e-05;
              result = -2.0019531249999998e-03 + result * t;
              result = 1.6547309027777776e-02 + result * t;
              result = -6.5104166666666671e-02 + result * t;
              result = 1.0850694444444445e-01 + result * t;
              result *= t;
              result = -1.4467592592592593e-01 + result * t;
              return innerDeriv * result;
            } else {
              t -= 3.0;
              double result = -4.6296296296296294e-05;
              result = 2.7777777777777778e-04 + result * t;
              result = -6.9444444444444447e-04 + result * t;
              result = 9.2592592592592596e-04 + result * t;
              result = -6.9444444444444447e-04 + result * t;
              result = 2.7777777777777778e-04 + result * t;
              result = -4.6296296296296294e-05 + result * t;
              return innerDeriv * result;
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
};

// default type-def (unsigned int for level and index)
typedef NakBsplineModifiedBasis<unsigned int, unsigned int> SNakBsplineModifiedBase;

}  // namespace base
}  // namespace sgpp