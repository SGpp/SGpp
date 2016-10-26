// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BSPLINE_NOT_A_KNOT_BOUNDARY_BASE_HPP
#define BSPLINE_NOT_A_KNOT_BOUNDARY_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace base {

/**
 * B-spline basis on Boundary grids with not-a-knot boundary conditions.
 */
template <class LT, class IT>
class BsplineNotAKnotBoundaryBasis: public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  BsplineNotAKnotBoundaryBasis() : bsplineBasis(BsplineBasis<LT, IT>()) {
  }

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit BsplineNotAKnotBoundaryBasis(size_t degree) :
    bsplineBasis(BsplineBasis<LT, IT>(degree)) {
    if (degree != 3) {
      throw std::out_of_range("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~BsplineNotAKnotBoundaryBasis() override {
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of boundary B-spline basis function
   */
  inline double eval(LT l, IT i, double x) override {
    if (l == 0) {
      if (i == 0) {
        return 1.0 - x;
      } else {
        return x;
      }
    } else if (l == 1) {
      return 4.0 * x * (1.0 - x);
    } else if (l == 2) {
      double t = 4.0 * x - static_cast<double>(i);

      if (i == 3) {
        t = -t;
      }

      if (t < 1.0) {
        t += 1.0;
        double y = 1.0/10.0;
        y = -9.0/20.0 + y * t;
        y = 3.0/10.0 + y * t;
        y = 3.0/5.0 + y * t;
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
      const double hInv = static_cast<double>(static_cast<IT>(1) << l);
      double t = x * hInv - static_cast<double>(i);

      if ((i == 1) || (i == hInv - 1)) {
        if (i == hInv - 1) {
          t = -t;
        }

        if (t < 1.0) {
          t += 1.0;
          double y = 1.0/8.0;
          y = -1.0/2.0 + y * t;
          y = 1.0/4.0 + y * t;
          y = 7.0/12.0 + y * t;
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
      } else if ((i == 3) || (i == hInv - 3)) {
        if (i == hInv - 3) {
          t = -t;
        }

        if (t < -1.0) {
          t += 3.0;
          double y = 1.0/24.0;
          y *= t;
          y *= t;
          y *= t;
          return y;
        } else if (t < 0.0) {
          t += 1.0;
          double y = -3.0/8.0;
          y = 1.0/4.0 + y * t;
          y = 1.0/2.0 + y * t;
          y = 1.0/3.0 + y * t;
          return y;
        } else if (t < 1.0) {
          double y = 11.0/24.0;
          y = -7.0/8.0 + y * t;
          y = -1.0/8.0 + y * t;
          y = 17.0/24.0 + y * t;
          return y;
        } else if (t < 2.0) {
          t -= 1.0;
          double y = -1.0/6.0;
          y = 1.0/2.0 + y * t;
          y = -1.0/2.0 + y * t;
          y = 1.0/6.0 + y * t;
          return y;
        } else {
          return 0.0;
        }
      } else {
        return bsplineBasis.uniformBSpline(
            t + static_cast<double>(bsplineBasis.getDegree() + 1) / 2.0,
            bsplineBasis.getDegree());
      }
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of derivative of boundary B-spline basis function
   */
  inline double evalDx(LT l, IT i, double x) {
    if (l == 0) {
      if (i == 0) {
        return -1.0;
      } else {
        return 1.0;
      }
    } else if (l == 1) {
      return 4.0 - 8.0 * x;
    } else if (l == 2) {
      double t = 4.0 * x - static_cast<double>(i);
      double inner_derivative = 4.0;

      if (i == 3) {
        t = -t;
        inner_derivative *= -1.0;
      }

      if (t < 1.0) {
        t += 1.0;
        double y = 3.0/10.0;
        y = -9.0/10.0 + y * t;
        y = 3.0/10.0 + y * t;
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
      double t = x * hInv - static_cast<double>(i);
      double inner_derivative = hInv;

      if ((i == 1) || (i == hInv - 1)) {
        if (i == hInv - 1) {
          t = -t;
          inner_derivative *= -1.0;
        }

        if (t < 1.0) {
          t += 1.0;
          double y = 3.0/8.0;
          y = -1.0 + y * t;
          y = 1.0/4.0 + y * t;
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
      } else if ((i == 3) || (i == hInv - 3)) {
        if (i == hInv - 3) {
          t = -t;
          inner_derivative *= -1.0;
        }

        if (t < -1.0) {
          t += 3.0;
          double y = 1.0/8.0;
          y *= t;
          y *= t;
          return inner_derivative * y;
        } else if (t < 0.0) {
          t += 1.0;
          double y = -9.0/8.0;
          y = 1.0/2.0 + y * t;
          y = 1.0/2.0 + y * t;
          return inner_derivative * y;
        } else if (t < 1.0) {
          double y = 11.0/8.0;
          y = -7.0/4.0 + y * t;
          y = -1.0/8.0 + y * t;
          return inner_derivative * y;
        } else if (t < 2.0) {
          t -= 1.0;
          double y = -1.0/2.0;
          y = 1.0 + y * t;
          y = -1.0/2.0 + y * t;
          return inner_derivative * y;
        } else {
          return 0.0;
        }
      } else {
        return inner_derivative * bsplineBasis.uniformBSplineDx(
            t + static_cast<double>(bsplineBasis.getDegree() + 1) / 2.0,
            bsplineBasis.getDegree());
      }
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of 2nd derivative of boundary B-spline basis function
   */
  inline double evalDxDx(LT l, IT i, double x) {
    if (l == 0) {
      return 0.0;
    } else if (l == 1) {
      return -8.0;
    } else if (l == 2) {
      double t = 4.0 * x - static_cast<double>(i);
      double inner_derivative = 4.0;

      if (i == 3) {
        t = -t;
        inner_derivative *= -1.0;
      }

      if (t < 1.0) {
        t += 1.0;
        double y = 3.0/5.0;
        y = -9.0/10.0 + y * t;
        return inner_derivative * inner_derivative * y;
      } else {
        t -= 1.0;
        double y = -3.0/20.0;
        y = 3.0/10.0 + y * t;
        return inner_derivative * inner_derivative * y;
      }
    } else {
      const double hInv = static_cast<double>(static_cast<IT>(1) << l);
      double t = x * hInv - static_cast<double>(i);
      double inner_derivative = hInv;

      if ((i == 1) || (i == hInv - 1)) {
        if (i == hInv - 1) {
          t = -t;
          inner_derivative *= -1.0;
        }

        if (t < 1.0) {
          t += 1.0;
          double y = 3.0/4.0;
          y = -1.0 + y * t;
          return inner_derivative * inner_derivative * y;
        } else if (t < 2.0) {
          t -= 1.0;
          double y = -1.0/2.0;
          y = 1.0/2.0 + y * t;
          return inner_derivative * inner_derivative * y;
        } else {
          return 0.0;
        }
      } else if ((i == 3) || (i == hInv - 3)) {
        if (i == hInv - 3) {
          t = -t;
          inner_derivative *= -1.0;
        }

        if (t < -1.0) {
          t += 3.0;
          double y = 1.0/4.0;
          y *= t;
          return inner_derivative * inner_derivative * y;
        } else if (t < 0.0) {
          t += 1.0;
          double y = -9.0/4.0;
          y = 1.0/2.0 + y * t;
          return inner_derivative * inner_derivative * y;
        } else if (t < 1.0) {
          double y = 11.0/4.0;
          y = -7.0/4.0 + y * t;
          return inner_derivative * inner_derivative * y;
        } else if (t < 2.0) {
          t -= 1.0;
          double y = -1.0;
          y = 1.0 + y * t;
          return inner_derivative * inner_derivative * y;
        } else {
          return 0.0;
        }
      } else {
        return inner_derivative * inner_derivative * bsplineBasis.uniformBSplineDxDx(
            t + static_cast<double>(bsplineBasis.getDegree() + 1) / 2.0,
            bsplineBasis.getDegree());
      }
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
};

// default type-def (unsigned int for level and index)
typedef BsplineNotAKnotBoundaryBasis<unsigned int, unsigned int> SBsplineNotAKnotBoundaryBase;

}  // namespace base
}  // namespace sgpp

#endif /* BSPLINE_NOT_A_KNOT_BOUNDARY_BASE_HPP */
