// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BSPLINE_CLENSHAW_CURTIS_BASE_HPP
#define BSPLINE_CLENSHAW_CURTIS_BASE_HPP

#include <cmath>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/tools/ClenshawCurtisTable.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace base {

/**
 * B-spline basis on Clenshaw-Curtis grids.
 */
template<class LT, class IT>
class BsplineClenshawCurtisBasis: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  BsplineClenshawCurtisBasis()
    : bsplineBasis(BsplineBasis<LT, IT>()) {
  }

  /**
   * Constructor.
   *
   * @param degree        B-spline degree, must be odd
   *                      (if it's even, degree - 1 is used)
   */
  BsplineClenshawCurtisBasis(size_t degree)
    : bsplineBasis(BsplineBasis<LT, IT>(degree)),
      xi(std::vector<float_t>(degree + 2, 0.0)) {
  }

  /**
   * Destructor.
   */
  virtual ~BsplineClenshawCurtisBasis() override {
  }

  /**
   * @param x     evaluation point
   * @param p     B-spline degree
   * @param k     index of B-spline in the knot sequence
   * @return      value of non-uniform B-spline
   *              with knots \f$\{\xi_k, ... \xi_{k+p+1}\}\f$
   */
  inline float_t nonUniformBSpline(float_t x, size_t p, size_t k) const {
    if (p == 0) {
      // characteristic function of [xi[k], xi[k+1])
      return (((x >= xi[k]) && (x < xi[k + 1])) ? 1.0 : 0.0);
    } else if ((x < xi[k]) || (x >= xi[k + p + 1])) {
      // out of support
      return 0.0;
    } else {
      // Cox-de-Boor recursion
      return (x - xi[k]) / (xi[k + p] - xi[k])
             * nonUniformBSpline(x, p - 1, k)
             + (1.0 - (x - xi[k + 1]) / (xi[k + p + 1] - xi[k + 1]))
             * nonUniformBSpline(x, p - 1, k + 1);
    }
  }

  /**
   * @param x     evaluation point
   * @param p     B-spline degree
   * @param k     index of B-spline in the knot sequence
   * @return      value of derivative of non-uniform B-spline
   *              with knots \f$\{\xi_k, ... \xi_{k+p+1}\}\f$
   */
  inline float_t nonUniformBSplineDx(
    float_t x, size_t p, size_t k) const {
    if (p == 0) {
      return 0.0;
    } else if ((x < xi[k]) || (x >= xi[k + p + 1])) {
      return 0.0;
    } else {
      const float_t pDbl = static_cast<float_t>(p);

      return pDbl / (xi[k + p] - xi[k]) * nonUniformBSpline(x, p - 1, k)
             - pDbl / (xi[k + p + 1] - xi[k + 1])
             * nonUniformBSpline(x, p - 1, k + 1);
    }
  }

  /**
   * @param x     evaluation point
   * @param p     B-spline degree
   * @param k     index of B-spline in the knot sequence
   * @return      value of 2nd derivative of non-uniform B-spline
   *              with knots \f$\{\xi_k, ... \xi_{k+p+1}\}\f$
   */
  inline float_t nonUniformBSplineDxDx(
    float_t x, size_t p, size_t k) const {
    if (p <= 1) {
      return 0.0;
    } else if ((x < xi[k]) || (x >= xi[k + p + 1])) {
      return 0.0;
    } else {
      const float_t pDbl = static_cast<float_t>(p);
      const float_t alphaKP = pDbl / (xi[k + p] - xi[k]);
      const float_t alphaKp1P = pDbl / (xi[k + p + 1] - xi[k + 1]);
      const float_t alphaKPm1 = (pDbl - 1.0) / (xi[k + p - 1] - xi[k]);
      const float_t alphaKp1Pm1 = (pDbl - 1.0) / (xi[k + p] - xi[k + 1]);
      const float_t alphaKp2Pm1 = (pDbl - 1.0) /
                                  (xi[k + p + 1] - xi[k + 2]);

      return alphaKP * alphaKPm1 * nonUniformBSpline(x, p - 2, k)
             - (alphaKP + alphaKp1P) * alphaKp1Pm1
             * nonUniformBSpline(x, p - 2, k + 1)
             + alphaKp1P * alphaKp2Pm1 *
             nonUniformBSpline(x, p - 2, k + 2);
    }
  }

  /**
   * @param l     level of the grid point
   * @param i     index of the grid point
   * @return      i-th Clenshaw-Curtis grid point with level l
   */
  inline float_t clenshawCurtisPoint(LT l, IT i) const {
    return ClenshawCurtisTable::getInstance().getPoint(l, i);
  }

  /**
   * Construct the (p+2) Clenshaw-Curtis knots of a
   * B-spline basis function and save them in xi.
   *
   * @param l     level of basis function
   * @param i     index of basis function
   */
  inline void constructKnots(LT l, IT i) {
    const IT hInv = static_cast<IT>(1) << l;
    const size_t& p = bsplineBasis.getDegree();

    xi[(p + 1) / 2] = ClenshawCurtisTable::getInstance().getPoint(l, i);

    if (i < (p + 1) / 2) {
      // grid point index is too far on the left
      // ==> extrapolate grid points linearly
      size_t a = (p + 1) / 2 - i;

      for (size_t j = a; j < (p + 1) / 2; j++) {
        xi[j] = ClenshawCurtisTable::getInstance().getPoint(l, static_cast<IT>(j - a));
      }

      float_t h = xi[a + 1] - xi[a];

      // equivalent to "for (int j = a-1; j >= 0; j--)"
      for (size_t j = a; j-- > 0;) {
        xi[j] = xi[j + 1] - h;
      }
    } else {
      // all grid points on the left can be calculated
      for (size_t j = 0; j < (p + 1) / 2; j++) {
        xi[j] = ClenshawCurtisTable::getInstance().getPoint(
                  l, static_cast<IT>(i - (p + 1) / 2 + j));
      }
    }

    if (i + (p + 1) / 2 > hInv) {
      // grid point index is too far on the right
      // ==> extrapolate grid points linearly
      size_t b = hInv + (p + 1) / 2 - i;

      for (size_t j = (p + 1) / 2 + 1; j <= b; j++) {
        xi[j] = ClenshawCurtisTable::getInstance().getPoint(
                  l, static_cast<IT>(i - (p + 1) / 2 + j));
      }

      float_t h = xi[b] - xi[b - 1];

      for (size_t j = b + 1; j < p + 2; j++) {
        xi[j] = xi[j - 1] + h;
      }
    } else {
      // all grid points on the right can be calculated
      for (size_t j = (p + 1) / 2 + 1; j < p + 2; j++) {
        xi[j] = ClenshawCurtisTable::getInstance().getPoint(
                  l, static_cast<IT>(i - (p + 1) / 2 + j));
      }
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of Clenshaw-Curtis B-spline basis function
   */
  inline virtual float_t eval(LT l, IT i, float_t x) override {
    if (l == 0) {
      return bsplineBasis.uniformBSpline(
               x - static_cast<float_t>(i)
               + static_cast<float_t>(bsplineBasis.getDegree() + 1) / 2.0,
               bsplineBasis.getDegree());
    } else {
      constructKnots(l, i);
      return nonUniformBSpline(x, bsplineBasis.getDegree(), 0);
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of derivative of Clenshaw-Curtis
   *              B-spline basis function
   */
  inline float_t evalDx(LT l, IT i, float_t x) {
    if (l == 0) {
      return bsplineBasis.uniformBSplineDx(
               x - static_cast<float_t>(i)
               + static_cast<float_t>(bsplineBasis.getDegree() + 1) / 2.0,
               bsplineBasis.getDegree());
    } else {
      constructKnots(l, i);
      return nonUniformBSplineDx(x, bsplineBasis.getDegree(), 0);
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of 2nd derivative of Clenshaw-Curtis
   *              B-spline basis function
   */
  inline float_t evalDxDx(LT l, IT i, float_t x) {
    if (l == 0) {
      return bsplineBasis.uniformBSplineDxDx(
               x - static_cast<float_t>(i)
               + static_cast<float_t>(bsplineBasis.getDegree() + 1) / 2.0,
               bsplineBasis.getDegree());
    } else {
      constructKnots(l, i);
      return nonUniformBSplineDxDx(x, bsplineBasis.getDegree(), 0);
    }
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const {
    return bsplineBasis.getDegree();
  }

 protected:
  /// B-spline basis for B-spline evaluation
  BsplineBasis<LT, IT> bsplineBasis;
  /// temporary helper vector of fixed size p+2 containing B-spline knots
  std::vector<float_t> xi;
};

// default type-def (unsigned int for level and index)
typedef BsplineClenshawCurtisBasis<unsigned int, unsigned int>
SBsplineClenshawCurtisBase;
}
}

#endif /* BSPLINE_CLENSHAW_CURTIS_BASE_HPP */
