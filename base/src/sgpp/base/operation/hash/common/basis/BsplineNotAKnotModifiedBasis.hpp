// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BSPLINE_NOT_A_KNOT_MODIFIED_BASE_HPP
#define BSPLINE_NOT_A_KNOT_MODIFIED_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineNotAKnotBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace base {

/**
 * Modified B-spline basis on Noboundary grids with not-a-knot boundary conditions.
 */
template <class LT, class IT>
class BsplineNotAKnotModifiedBasis: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  BsplineNotAKnotModifiedBasis() :
    bsplineNotAKnotBoundaryBasis(BsplineNotAKnotBoundaryBasis<LT, IT>()) {
  }

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit BsplineNotAKnotModifiedBasis(size_t degree) :
    bsplineNotAKnotBoundaryBasis(BsplineNotAKnotBoundaryBasis<LT, IT>(degree)) {
    if (degree != 3) {
      throw std::out_of_range("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~BsplineNotAKnotModifiedBasis() override {
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of modified B-spline basis function
   */
  inline double eval(LT l, IT i, double x) override {
    if (l == 1) {
      return 1.0;
    } else if (l == 2) {
      double t = 4.0 * x - static_cast<double>(i);

      if (i == 3) {
        t = -t;
      }

      if (t < 1.0) {
        t += 1.0;
        double y = 1.0/40.0;
        y *= t;
        y = -3.0/5.0 + y * t;
        y = 6.0/5.0 + y * t;
        return y;
      } else {
        t -= 1.0;
        double y = -1.0/40.0;
        y = 3.0/20.0 + y * t;
        y = -3.0/10.0 + y * t;
        y = 1.0/5.0 + y * t;
        return y;
      }
    } else {
      double hInv = static_cast<double>(static_cast<IT>(1) << l);

      if ((i == 1) || (i == hInv - 1)) {
        double t = x * hInv - static_cast<double>(i);

        if (i == hInv - 1) {
          t = -t;
        }

        if (t < 1.0) {
          t += 1.0;
          double y = 1.0/24.0;
          y *= t;
          y = -3.0/4.0 + y * t;
          y = 5.0/4.0 + y * t;
          return y;
        } else if (t < 2.0) {
          t -= 1.0;
          double y = -1.0/12.0;
          y = 1.0/4.0 + y * t;
          y = -1.0/4.0 + y * t;
          y = 1.0/12.0 + y * t;
          return y;
        } else {
          return 0.0;
        }
      } else {
        return bsplineNotAKnotBoundaryBasis.eval(l, i, x);
      }
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of derivative of modified
   *              B-spline basis function
   */
  inline double evalDx(LT l, IT i, double x) {
    if (l == 1) {
      return 0.0;
    } else if (l == 2) {
      double t = 4.0 * x - static_cast<double>(i);
      double inner_derivative = 4.0;

      if (i == 3) {
        t = -t;
        inner_derivative *= -1.0;
      }

      if (t < 1.0) {
        t += 1.0;
        double y = 3.0/40.0;
        y *= t;
        y = -3.0/5.0 + y * t;
        return inner_derivative * y;
      } else {
        t -= 1.0;
        double y = -3.0/40.0;
        y = 3.0/10.0 + y * t;
        y = -3.0/10.0 + y * t;
        return inner_derivative * y;
      }
    } else {
      const double hInv = static_cast<double>(static_cast<IT>(1) << l);

      if ((i == 1) || (i == hInv - 1)) {
        double t = x * hInv - static_cast<double>(i);
        double inner_derivative = hInv;

        if (i == hInv - 1) {
          t = -t;
          inner_derivative *= -1.0;
        }

        if (t < 1.0) {
          t += 1.0;
          double y = 1.0/8.0;
          y *= t;
          y = -3.0/4.0 + y * t;
          return inner_derivative * y;
        } else if (t < 2.0) {
          t -= 1.0;
          double y = -1.0/4.0;
          y = 1.0/2.0 + y * t;
          y = -1.0/4.0 + y * t;
          return inner_derivative * y;
        } else {
          return 0.0;
        }
      } else {
        return bsplineNotAKnotBoundaryBasis.evalDx(l, i, x);
      }
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of 2nd derivative of modified
   *              B-spline basis function
   */
  inline double evalDxDx(LT l, IT i, double x) {
    if (l == 1) {
      return 0.0;
    } else if (l == 2) {
      double t = 4.0 * x - static_cast<double>(i);
      double inner_derivative = 4.0;

      if (i == 3) {
        t = -t;
        inner_derivative *= -1.0;
      }

      if (t < 1.0) {
        t += 1.0;
        double y = 3.0/20.0;
        y *= t;
        return inner_derivative * inner_derivative * y;
      } else {
        t -= 1.0;
        double y = -3.0/20.0;
        y = 3.0/10.0 + y * t;
        return inner_derivative * inner_derivative * y;
      }
    } else {
      const double hInv = static_cast<double>(static_cast<IT>(1) << l);

      if ((i == 1) || (i == hInv - 1)) {
        double t = x * hInv - static_cast<double>(i);
        double inner_derivative = hInv;

        if (i == hInv - 1) {
          t = -t;
          inner_derivative *= -1.0;
        }

        if (t < 1.0) {
          t += 1.0;
          double y = 1.0/4.0;
          y *= t;
          return inner_derivative * inner_derivative * y;
        } else if (t < 2.0) {
          t -= 1.0;
          double y = -1.0/2.0;
          y = 1.0/2.0 + y * t;
          return inner_derivative * inner_derivative * y;
        } else {
          return 0.0;
        }
      } else {
        return bsplineNotAKnotBoundaryBasis.evalDxDx(l, i, x);
      }
    }
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const {
    return bsplineNotAKnotBoundaryBasis.getDegree();
  }

 protected:
  /// B-spline not-a-knot boundary basis for B-spline evaluation
  BsplineNotAKnotBoundaryBasis<LT, IT> bsplineNotAKnotBoundaryBasis;
};

// default type-def (unsigned int for level and index)
typedef BsplineNotAKnotModifiedBasis<unsigned int, unsigned int>
SBsplineNotAKnotModifiedBase;

}  // namespace base
}  // namespace sgpp

#endif /* BSPLINE_NOT_A_KNOT_MODIFIED_BASE_HPP */
