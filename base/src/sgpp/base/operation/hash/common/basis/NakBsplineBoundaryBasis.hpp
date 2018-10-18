// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>

namespace sgpp {
namespace base {

/**
 * Hierarchical Not-a-knot B-spline basis.
 *
 * This basis has the canonical choice of a constant basis function on level 0 in contrast to
 * NakBsplineBoundaryCombigridBasis which does otherwise to become analogous to the combigrid not a
 * knot B spline boundary basis.
 */

template <class LT, class IT>
class NakBsplineBoundaryBasis : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakBsplineBoundaryBasis() : nakBsplineBasis(NakBsplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineBoundaryBasis(size_t degree)
      : nakBsplineBasis(NakBsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 7) {
      throw std::runtime_error("NakBsplineBoundaryBasis: Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakBsplineBoundaryBasis() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double eval(LT l, IT i, double x) override {
    switch (getDegree()) {
      case 1:
        return nakBsplineBasis.eval(l, i, x);

      // degree 3: global polynomials in x on Level 0 and 1, nak Bsplines from Level 3 on
      case 3:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            // return 1;
            return 1 - x;
          } else {
            // l = 0, i = 1
            return x;
          }
        } else if (l == 1) {
          if (i == 0) {
            // l = 1, i = 0
            return -2 * (1 - x) * (x - 0.5);
          } else if (i == 1) {
            // l = 1, i = 1
            // return -4 * (x - 1) * x;
            return 4 * x - 4 * x * x;
          } else {
            // l = 1, i = 2
            return 2 * (x - 0.5) * x;
          }

        } else {
          return nakBsplineBasis.eval(l, i, x);
        }

      // degree 5: Levels 0,1 and 2 polynomials, nak Bsplines from Level 3 on
      case 5:

        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return 1;
          } else {
            // l = 0, i = 1
            return x;
          }
        } else if (l == 1) {
          if (i == 1) {
            // l = 1, i = 1
            return -4 * (x - 1) * x;
          }
        } else if (l == 2) {
          if (i == 1) {
            // l = 2, i = 1 : cubic polynomial, 0 in 0,0.5,0.75 and 1 in 0.25
            return 32 * x * (x - 0.5) * (x - 0.75);
          } else if (i == 3) {
            // l = 2, i = 3 : quartic polynomial, 0 in 0,0.25,0.5,1 and 1 in 0.75
            return x * x * x * x;  // x * (x - 0.25) * (x - 0.5) * (x - 1) * (-128.0 / 3.0);
          }
        } else {
          return nakBsplineBasis.eval(l, i, x);
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
   * @param x     evaluation point
   * @return      value of derivative of basis function
   */
  inline double evalDx(LT l, IT i, double x) {
    switch (getDegree()) {
      case 1:
        return nakBsplineBasis.evalDx(l, i, x);
      case 3:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return -1;
          } else {
            // l = 0, i = 1
            return 1;
          }
        } else if (l == 1) {
          if (i == 0) {
            // l = 1, i = 0
            return 4 * (x - 0.75);
          } else if (i == 1) {
            // l = 1, i = 1
            // return -4 * (x - 1) * x;
            return 4 - 8 * x;
          } else {
            // l = 1, i = 2
            return 4 * (x - 0.25);
          }

        } else {
          return nakBsplineBasis.evalDx(l, i, x);
        }
      case 5:
        if (l == 0) {
          if (i == 0) {
            // l = 0, i = 0
            return 0;
          } else {
            // l = 0, i = 1
            return 1;
          }
        } else if (l == 1) {
          if (i == 1) {
            // l = 1, i = 1
            return 4 - 8 * x;
          }
        } else if (l == 2) {
          if (i == 1) {
            // l = 2, i = 1
            return 96 * (x * x - 0.83333333333333333 * x + 0.125);
          } else if (i == 3) {
            // l = 2, i = 3
            return 4 * x * x * x;
          }
        } else {
          return nakBsplineBasis.eval(l, i, x);
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
  inline double getWeightedIntegral(LT l, IT i, sgpp::combigrid::SingleFunction weightfunction,
                                    double lbound, double rbound, size_t numAdditionalPoints = 0,
                                    size_t incrementQuadraturePoints = 10) {
    size_t degree = getDegree();
    if ((degree != 1) && (degree != 3) && (degree != 5)) {
      throw std::runtime_error(
          "OperationMatrixLTwoDotNakBsplineBoundary: only B spline degrees 1, 3 and 5 are "
          "supported.");
    }
    size_t quadOrder = degree + 1 + numAdditionalPoints;
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

    double temp_res = integrateWeightedBspline(l, i, start, stop, offset, scaling, quadCoordinates,
                                               quadWeights, weightfunction);
    double width = rbound - lbound;
    temp_res *= width;

    double tol = 1e-14;
    double err = 1e14;
    while (err > tol) {
      numAdditionalPoints += incrementQuadraturePoints;
      quadOrder = degree + 1 + numAdditionalPoints;
      if (quadOrder > 480) {
        break;
      }
      gauss.getLevelPointsAndWeightsNormalized(quadOrder, quadCoordinates, quadWeights);
      double finer_temp_res = integrateWeightedBspline(
          l, i, start, stop, offset, scaling, quadCoordinates, quadWeights, weightfunction);
      finer_temp_res *= width;
      err = fabs(temp_res - finer_temp_res);
      temp_res = finer_temp_res;
    }

    double integral = temp_res * scaling;
    return integral;
  }
  //#endif

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const { return nakBsplineBasis.getDegree(); }

 protected:
  /// not a knot B-spline basis
  NakBsplineBasis<LT, IT> nakBsplineBasis;

 private:
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
  double integrateWeightedBspline(LT l, IT i, size_t start, size_t stop, double offset,
                                  double scaling, base::DataVector quadCoordinates,
                                  base::DataVector quadWeights,
                                  sgpp::combigrid::SingleFunction weightfunction) {
    double temp_res = 0.0;
    for (size_t n = start; n <= stop; n++) {
      for (size_t c = 0; c < quadCoordinates.getSize(); c++) {
        // transform  the quadrature points to the segment on which the Bspline is
        // evaluated
        const double x = offset + scaling * (quadCoordinates[c] + static_cast<double>(n));
        temp_res += quadWeights[c] * this->eval(l, i, x) * weightfunction(x);
      }
    }
    return temp_res;
  }
  //#endif
};

// default type-def (unsigned int for level and index)
typedef NakBsplineBoundaryBasis<unsigned int, unsigned int> SNakBsplineBoundaryBase;

}  // namespace base
}  // namespace sgpp
