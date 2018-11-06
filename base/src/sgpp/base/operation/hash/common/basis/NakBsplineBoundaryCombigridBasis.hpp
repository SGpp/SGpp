// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
//#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>

namespace sgpp {
namespace base {

/**
 * Hierarchical Not-a-knot B-spline basis. This basis is designed to represent Bspline boundary
 * interpolants from the combigrid module.
 * Therefore it has the following unusual choice of basis functions:
 * linear and quadratic terms on level 0
 * constant term on level 1
 *
 * This means that this basis of level 0 cannot represent constant functions! It is therefore not
 * suitabe for any application that cannot guarantee the existence of level 1 in every dimension
 *
 * A combigrid interpolant which shall be interpolated with this class needs level (1,...,1)
 * otherwise the boundary points cannot match
 *
 * If the interpolant has not been created by the combigrid module use the NakBsplineBoundaryBasis
 * instead
 */

template <class LT, class IT>
class NakBsplineBoundaryCombigridBasis : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakBsplineBoundaryCombigridBasis() : bsplineBasis(BsplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineBoundaryCombigridBasis(size_t degree)
      : bsplineBasis(BsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("NakBsplineBoundaryCombigridBasis: Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakBsplineBoundaryCombigridBasis() override {}

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

      // degree 3: global polynomials in x on Level 0 and 1, nak Bsplines from Level 3 on
      case 3:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            //            return 1 - x * x;
            return (2 * x - 3) * x + 1;
          } else {
            // l = 0, i = 1
            //            return x;
            return (2 * x - 1) * x;
          }
        } else if (l == 1) {
          if (i == 0) {
            // l = 1, i = 0
            return 0.5 * x * x - 1.5 * x + 1.0;
          } else if (i == 1) {
            // l = 1, i = 1
            return 1;
          } else {
            // l = 1, i = 2
            return 0.5 * x * x + 1.5 * x + 1.0;
          }
        } else if ((i > 3) && (i < hInv - 3)) {
          // l >= 4, 3 < i < 2^l - 3
          return bsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if (i == 0) {
            // l >= 2, i = 0
            if ((t < 0.0) || (t > 2.0)) {
              return 0.0;
            } else {
              double result = -4.1666666666666664e-02;
              result = 2.5000000000000000e-01 + result * t;
              result = -5.0000000000000000e-01 + result * t;
              result = 3.3333333333333331e-01 + result * t;
              return result;
            }
          } else if ((l == 2) && (i == 1)) {
            // l = 2, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 1.0 / 10.0;
              result = -9.0 / 20.0 + result * t;
              result = 3.0 / 10.0 + result * t;
              result = 3.0 / 5.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0 / 40.0;
              result = 3.0 / 20.0 + result * t;
              result = -3.0 / 10.0 + result * t;
              result = 1.0 / 5.0 + result * t;
              return result;
            }

          } else if (l == 2) {
            // l = 2, i = 2
            if ((t < -2.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 2.0;
              double result = -8.3333333333333329e-02;
              result = 2.0000000000000001e-01 + result * t;
              result = 2.0000000000000001e-01 + result * t;
              result = 6.6666666666666666e-02 + result * t;
              return result;
            } else {
              double result = 8.3333333333333329e-02;
              result = -2.9999999999999999e-01 + result * t;
              result *= t;
              result = 5.9999999999999998e-01 + result * t;
              return result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 1.0) {
              t += 1.0;
              double result = 1.0 / 8.0;
              result = -1.0 / 2.0 + result * t;
              result = 1.0 / 4.0 + result * t;
              result = 7.0 / 12.0 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.0 / 12.0;
              result = 1.0 / 4.0 + result * t;
              result = -1.0 / 4.0 + result * t;
              result = 1.0 / 12.0 + result * t;
              return result;
            }
          } else if (i == 2) {
            // l >= 3, i = 2
            if ((t < -2.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 2.0;
              double result = -1.2500000000000000e-01;
              result = 2.5000000000000000e-01 + result * t;
              result = 2.5000000000000000e-01 + result * t;
              result = 8.3333333333333329e-02 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 2.9166666666666669e-01;
              result = -5.0000000000000000e-01 + result * t;
              result = -2.5000000000000000e-01 + result * t;
              result = 5.8333333333333337e-01 + result * t;
              return result;
            } else {
              t -= 1.0;
              double result = -1.2500000000000000e-01;
              result = 3.7500000000000000e-01 + result * t;
              result = -3.7500000000000000e-01 + result * t;
              result = 1.2500000000000000e-01 + result * t;
              return result;
            }
          } else {
            // l >= 3, i = 3
            if ((t < -3.0) || (t > 2.0)) {
              return 0.0;
            } else if (t < -1.0) {
              t += 3.0;
              double result = 1.0 / 24.0;
              result *= t;
              result *= t;
              result *= t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = -3.0 / 8.0;
              result = 1.0 / 4.0 + result * t;
              result = 1.0 / 2.0 + result * t;
              result = 1.0 / 3.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = 11.0 / 24.0;
              result = -7.0 / 8.0 + result * t;
              result = -1.0 / 8.0 + result * t;
              result = 17.0 / 24.0 + result * t;
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

      // degree 5: Levels 0,1 and 2 polynomials, nak Bsplines from Level 3 on
      case 5:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return 1 - x * x;
          } else {
            // l = 0, i = 1
            return x;
          }
        } else if (l == 1) {
          if (i == 1) {
            // l = 1, i = 1
            return 1;
          }
        } else if (l == 2) {
          if (i == 1) {
            // l = 2, i = 1 : cubic polynomial, 0 in 0,0.5,0.75 and 1 in 0.25
            return 32 * x * (x - 0.5) * (x - 0.75);
          } else if (i == 3) {
            // l = 2, i = 3 : quartic polynomial, 0 in 0,0.25,0.5,1 and 1 in 0.75
            return x * x * x * x;  // x * (x - 0.25) * (x - 0.5) * (x - 1) * (-128.0 / 3.0);
          }
        } else if ((i > 5) && (i < hInv - 5)) {
          // l >= 4, 5 < i < 2^l - 5
          return bsplineBasis.eval(l, i, x);
        } else {
          if (i > hInv / 2) {
            i = hInv - i;
            t *= -1.0;
          }

          if ((l == 3) && (i == 3)) {
            // l = 3, i = 3
            if ((t < -3.0) || (t > 5.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 107.0 / 30240.0;
              result = -17.0 / 756.0 + result * t;
              result = 1.0 / 378.0 + result * t;
              result = 37.0 / 378.0 + result * t;
              result = 109.0 / 756.0 + result * t;
              result = 253.0 / 3780.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -397.0 / 30240.0;
              result = 185.0 / 6048.0 + result * t;
              result = 155.0 / 3024.0 + result * t;
              result = -415.0 / 3024.0 + result * t;
              result = -1165.0 / 6048.0 + result * t;
              result = 2965.0 / 6048.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 233.0 / 30240.0;
              result = -53.0 / 1512.0 + result * t;
              result = 8.0 / 189.0 + result * t;
              result = 13.0 / 189.0 + result * t;
              result = -97.0 / 378.0 + result * t;
              result = 433.0 / 1890.0 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.0 / 4320.0;
              result = 1.0 / 288.0 + result * t;
              result = -1.0 / 48.0 + result * t;
              result = 1.0 / 16.0 + result * t;
              result = -3.0 / 32.0 + result * t;
              result = 9.0 / 160.0 + result * t;
              return result;
            }
          } else if (i == 1) {
            // l >= 3, i = 1
            if ((t < -1.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 2.0) {
              t += 1.0;
              double result = 1.0 / 504.0;
              result = -1.0 / 42.0 + result * t;
              result = 2.0 / 21.0 + result * t;
              result = -2.0 / 21.0 + result * t;
              result = -5.0 / 21.0 + result * t;
              result = 47.0 / 105.0 + result * t;
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
          } else if (i == 3) {
            // l >= 4, i = 3
            if ((t < -3.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < 0.0) {
              t += 3.0;
              double result = 1.0 / 252.0;
              result = -1.0 / 42.0 + result * t;
              result *= t;
              result = 2.0 / 21.0 + result * t;
              result = 1.0 / 7.0 + result * t;
              result = 1.0 / 15.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -23.0 / 1260.0;
              result = 1.0 / 28.0 + result * t;
              result = 1.0 / 14.0 + result * t;
              result = -5.0 / 42.0 + result * t;
              result = -1.0 / 4.0 + result * t;
              result = 163.0 / 420.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 19.0 / 1260.0;
              result = -1.0 / 18.0 + result * t;
              result = 2.0 / 63.0 + result * t;
              result = 8.0 / 63.0 + result * t;
              result = -2.0 / 9.0 + result * t;
              result = 34.0 / 315.0 + result * t;
              return result;
            } else {
              t -= 2.0;
              double result = -1.0 / 252.0;
              result = 5.0 / 252.0 + result * t;
              result = -5.0 / 126.0 + result * t;
              result = 5.0 / 126.0 + result * t;
              result = -5.0 / 252.0 + result * t;
              result = 1.0 / 252.0 + result * t;
              return result;
            }
          } else {
            // l >= 4, i = 5
            if ((t < -5.0) || (t > 3.0)) {
              return 0.0;
            } else if (t < -2.0) {
              t += 5.0;
              double result = 1.0 / 2520.0;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              result *= t;
              return result;
            } else if (t < -1.0) {
              t += 2.0;
              double result = -11.0 / 504.0;
              result = 1.0 / 168.0 + result * t;
              result = 1.0 / 28.0 + result * t;
              result = 3.0 / 28.0 + result * t;
              result = 9.0 / 56.0 + result * t;
              result = 27.0 / 280.0 + result * t;
              return result;
            } else if (t < 0.0) {
              t += 1.0;
              double result = 31.0 / 504.0;
              result = -13.0 / 126.0 + result * t;
              result = -10.0 / 63.0 + result * t;
              result = 2.0 / 63.0 + result * t;
              result = 25.0 / 63.0 + result * t;
              result = 121.0 / 315.0 + result * t;
              return result;
            } else if (t < 1.0) {
              double result = -181.0 / 2520.0;
              result = 103.0 / 504.0 + result * t;
              result = 11.0 / 252.0 + result * t;
              result = -113.0 / 252.0 + result * t;
              result = -61.0 / 504.0 + result * t;
              result = 1543.0 / 2520.0 + result * t;
              return result;
            } else if (t < 2.0) {
              t -= 1.0;
              double result = 11.0 / 280.0;
              result = -13.0 / 84.0 + result * t;
              result = 1.0 / 7.0 + result * t;
              result = 4.0 / 21.0 + result * t;
              result = -3.0 / 7.0 + result * t;
              result = 23.0 / 105.0 + result * t;
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
        return 0.0;
        break;

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
          "OperationMatrixLTwoDotNakBsplineBoundary: only B spline degrees 1, 3 and 5 are "
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
   * @param numAdditionalPoints	number of additional points for the integration of the
   * weighted spline
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

  double evalDx(LT level, IT index, double x) override {
    std::cerr << "NakBsplineBoundaryCombigridBasis: evalDx not implemented" << std::endl;
    return -1;
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override { return bsplineBasis.getDegree(); }

 protected:
  /// B-spline basis for B-spline evaluation
  BsplineBasis<LT, IT> bsplineBasis;

 private:
  double integrateBspline(LT l, IT i, size_t start, size_t stop, double offset, double scaling,
                          base::DataVector quadCoordinates, base::DataVector quadWeights) {
    double temp_res = 0.0;
    for (size_t n = start; n <= stop; n++) {
      for (size_t c = 0; c < quadCoordinates.getSize(); c++) {
        // transform  the quadrature points to the segment on which the Bspline is
        // evaluated
        const double x = offset + scaling * (quadCoordinates[c] + static_cast<double>(n));
        if (this->eval(l, i, x) == 0) {
        }
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
typedef NakBsplineBoundaryCombigridBasis<unsigned int, unsigned int>
    SNakBsplineBoundaryCombigridBase;

}  // namespace base
}  // namespace sgpp
