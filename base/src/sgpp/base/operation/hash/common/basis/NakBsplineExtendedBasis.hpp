// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/tools/Distribution.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>

namespace sgpp {
namespace base {

/**
 * Extended Not-a-knot B-spline basis.
 *
 * This basis is a nonboundary version for the not-a-knot Bspline basis which still can represent
 * 1,x,...x^degree exact on [0,1].
 * This is achieved by adding the B-splines belonging to the left(0) and right boundary point (1) to
 * (degree+1) inner B-splines with extension coefficients derived from SLE for the representation of
 * 1,x,...x^degree.
 *
 * Extension of the not hierarchical basis:
 * For each of the outer B-splines b_0 and b_{2^l} a set of (degree+1) closest inner B-splines
 * I(j) is chosen. I.e. I(0) = {1,..,degree+1}, I(2^l) = {2^l-(degree+1),...,2^l-1}.
 * Now the functions 1,x,...,x^degree are interpolated with the not a knot B-splines basis resulting
 * in coefficients a^1_i,...,a^{x^degree}_i, i=0,...,2^l
 * Then the coefficients a^1_j,...,a^{x^degree}_j are represented as a linear combination of the
 * coefficients of I(j) for j in {0,2^l}.
 * The coefficients e_ij of these linear combinations are the extension coefficients
 *
 * We build a hierarchical basis from this by choosing Lagrange polynomials on levels 1 (and 2 for
 * degree 3 and 5) and then simply using the extended basis funtions from the nonhierarchical basis.
 *
 * (Note that one could also use the Lagrange polynomials in the sets I(0),I(8). This leads to
 * another spanned space. Also on can only extend once on level 3, then 1,x,...,x^degree are already
 * part of the spanned space and not extend on levels >=3 again leading to another spanned space.
 * ToDo (rehmemk) It might be interesting to compare the various possibilities of extended not a
 * knot B-spline bases
 * (Note: coefficients were calculated with ExtensionBspline.m))
 *
 */
template <class LT, class IT>
class NakBsplineExtendedBasis : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  NakBsplineExtendedBasis() : notAKnotBsplineBasis(NakBsplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit NakBsplineExtendedBasis(size_t degree)
      : notAKnotBsplineBasis(NakBsplineBasis<LT, IT>(degree)) {
    if (getDegree() > 5) {
      throw std::runtime_error("Unsupported B-spline degree.");
    }
  }

  /**
   * Destructor.
   */
  ~NakBsplineExtendedBasis() override {}

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;

    switch (getDegree()) {
        // Extension coefficients are left =[2, -1], right =[-1,2]
        // Only the ones corresponding to B-splines which are in the basis are used, i.e. indices
        // 1 and end-1

      case 1:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else {
          if (i == 1) {
            return notAKnotBsplineBasis.eval(l, i, x) + 2 * notAKnotBsplineBasis.eval(l, 0, x);
          } else if (i == hInv - 1) {
            return notAKnotBsplineBasis.eval(l, i, x) + 2 * notAKnotBsplineBasis.eval(l, hInv, x);
          } else {
            return notAKnotBsplineBasis.eval(l, i, x);
          }
        }

        // Constant on level 1
        // Lagrange Polynomials on level 2
        // Extended B-splines from level 3 on.

        // Extension coefficients are left =[5, -10, 10, -4], right =[-4,10,-10,5]
        // Only the ones corresponding to B-splines which are in the basis are used, i.e. indices
        // 1,3 and end-3,end-1
      case 3:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else if (l == 2) {
          if (i == 1) {
            // l = 2, i = 1
            return 8 * (x - 0.5) * (x - 0.75);
          } else {
            // l = 2, i = 3
            return 8 * (x - 0.25) * (x - 0.5);
          }

        } else {
          if (i == 1) {
            return notAKnotBsplineBasis.eval(l, i, x) + 5 * notAKnotBsplineBasis.eval(l, 0, x);
          } else if (i == 3) {
            return notAKnotBsplineBasis.eval(l, i, x) + 10 * notAKnotBsplineBasis.eval(l, 0, x);
          } else if (i == hInv - 3) {
            return notAKnotBsplineBasis.eval(l, i, x) + 10 * notAKnotBsplineBasis.eval(l, hInv, x);
          } else if (i == hInv - 1) {
            return notAKnotBsplineBasis.eval(l, i, x) + 5 * notAKnotBsplineBasis.eval(l, hInv, x);
          } else {
            return notAKnotBsplineBasis.eval(l, i, x);
          }
        }

        // Constant on level 1
        // Lagrange Polynomials on level 2
        // Extended B-splines from level 3 on.

        // Extension coefficients are left =[8,-28,42,-35,20,-6], right =[-6,20,-35,42,-28,8] on
        // level 3 and left = [8,-28,56,-70,56,-21], right = [-21,56,-70,56,-28,8] on level >= 4
        // Only the ones corresponding to B-splines which are in the basis are used, i.e. indices
        // 1,3,5 and end-5,end-3,end-1
      case 5:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else if (l == 2) {
          if (i == 1) {
            // l = 2, i = 1
            return 8 * (x - 0.5) * (x - 0.75);
          } else {
            // l = 2, i = 3
            return 8 * (x - 0.25) * (x - 0.5);
          }

        } else if (l == 3) {
          if (i == 1) {
            return notAKnotBsplineBasis.eval(l, i, x) + 8 * notAKnotBsplineBasis.eval(l, 0, x);
          } else if (i == 3) {
            return notAKnotBsplineBasis.eval(l, i, x) + 42 * notAKnotBsplineBasis.eval(l, 0, x) +
                   20 * notAKnotBsplineBasis.eval(l, hInv, x);
          } else if (i == 5) {
            return notAKnotBsplineBasis.eval(l, i, x) + 20 * notAKnotBsplineBasis.eval(l, 0, x) +
                   42 * notAKnotBsplineBasis.eval(l, hInv, x);
          } else if (i == 7) {
            return notAKnotBsplineBasis.eval(l, i, x) + 8 * notAKnotBsplineBasis.eval(l, hInv, x);
          } else {
            return notAKnotBsplineBasis.eval(l, i, x);
          }
        } else {
          if (i == 1) {
            return notAKnotBsplineBasis.eval(l, i, x) + 8 * notAKnotBsplineBasis.eval(l, 0, x);
          } else if (i == 3) {
            return notAKnotBsplineBasis.eval(l, i, x) + 56 * notAKnotBsplineBasis.eval(l, 0, x);
          } else if (i == 5) {
            return notAKnotBsplineBasis.eval(l, i, x) + 56 * notAKnotBsplineBasis.eval(l, 0, x);
          } else if (i == hInv - 5) {
            return notAKnotBsplineBasis.eval(l, i, x) + 56 * notAKnotBsplineBasis.eval(l, hInv, x);
          } else if (i == hInv - 3) {
            return notAKnotBsplineBasis.eval(l, i, x) + 56 * notAKnotBsplineBasis.eval(l, hInv, x);
          } else if (i == hInv - 1) {
            return notAKnotBsplineBasis.eval(l, i, x) + 8 * notAKnotBsplineBasis.eval(l, hInv, x);
          } else {
            return notAKnotBsplineBasis.eval(l, i, x);
          }
        }

      default:
        return 0.0;
    }
  }  // namespace base

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of basis function
   */
  inline double evalDx(LT l, IT i, double x) override {
    const IT hInv = static_cast<IT>(1) << l;
    switch (getDegree()) {
        // Extension coefficients are left =[2, -1], right =[-1,2]
        // Only the ones corresponding to B-splines which are in the basis are used, i.e. indices
        // 1,end-1

      case 1:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else {
          if (i == 1) {
            return notAKnotBsplineBasis.evalDx(l, i, x) + 2 * notAKnotBsplineBasis.evalDx(l, 0, x);
          } else if (i == hInv - 1) {
            return notAKnotBsplineBasis.evalDx(l, i, x) +
                   2 * notAKnotBsplineBasis.evalDx(l, hInv, x);
          } else {
            return notAKnotBsplineBasis.evalDx(l, i, x);
          }
        }

        // Constant on level 1
        // Lagrange Polynomials on level 2
        // Extended B-splines from level 3 on.

        // Extension coefficients are left =[5, -10, 10, -4], right =[-4,10,-10,5]
        // Only the ones corresponding to B-splines which are in the basis are used, i.e. indices
        // 1,3,end-3,end-1
      case 3:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if (l == 2) {
          if (i == 1) {
            // l = 2, i = 1
            return 16 * (x - 0.625);
          } else {
            // l = 2, i = 3
            return 16 * (x - 0.375);
          }

        } else {
          if (i == 1) {
            return notAKnotBsplineBasis.evalDx(l, i, x) + 5 * notAKnotBsplineBasis.evalDx(l, 0, x);
          } else if (i == 3) {
            return notAKnotBsplineBasis.evalDx(l, i, x) + 10 * notAKnotBsplineBasis.evalDx(l, 0, x);
          } else if (i == hInv - 3) {
            return notAKnotBsplineBasis.evalDx(l, i, x) +
                   10 * notAKnotBsplineBasis.evalDx(l, hInv, x);
          } else if (i == hInv - 1) {
            return notAKnotBsplineBasis.evalDx(l, i, x) +
                   5 * notAKnotBsplineBasis.evalDx(l, hInv, x);
          } else {
            return notAKnotBsplineBasis.evalDx(l, i, x);
          }
        }

        // Constant on level 1
        // Lagrange Polynomials on level 2
        // Extended B-splines from level 3 on.

        // Extension coefficients are left =[8,-28,42,-35,20,-6], right =[-6,20,-35,42,-28,8] on
        // level 3 and left = [8,-28,56,-70,56,-21], right = [-21,56,-70,56,-28,8] on level >= 4
        // Only the ones corresponding to B-splines which are in the basis are used, i.e. indices
        // 1,3,5,end-5,end-3,end-1
      case 5:
        if (l == 1) {
          // l = 1, i = 1
          return 0.0;
        } else if (l == 2) {
          if (i == 1) {
            // l = 2, i = 1
            return 16 * (x - 0.625);
          } else {
            // l = 2, i = 3
            return 16 * (x - 0.375);
          }

        } else if (l == 3) {
          if (i == 1) {
            return notAKnotBsplineBasis.evalDx(l, i, x) + 8 * notAKnotBsplineBasis.evalDx(l, 0, x);
          } else if (i == 3) {
            return notAKnotBsplineBasis.evalDx(l, i, x) +
                   42 * notAKnotBsplineBasis.evalDx(l, 0, x) +
                   20 * notAKnotBsplineBasis.evalDx(l, hInv, x);
          } else if (i == 5) {
            return notAKnotBsplineBasis.evalDx(l, i, x) +
                   20 * notAKnotBsplineBasis.evalDx(l, 0, x) +
                   42 * notAKnotBsplineBasis.evalDx(l, hInv, x);
          } else if (i == 7) {
            return notAKnotBsplineBasis.evalDx(l, i, x) +
                   8 * notAKnotBsplineBasis.evalDx(l, hInv, x);
          } else {
            return notAKnotBsplineBasis.evalDx(l, i, x);
          }
        } else {
          if (i == 1) {
            return notAKnotBsplineBasis.evalDx(l, i, x) + 8 * notAKnotBsplineBasis.evalDx(l, 0, x);
          } else if (i == 3) {
            return notAKnotBsplineBasis.evalDx(l, i, x) + 56 * notAKnotBsplineBasis.evalDx(l, 0, x);
          } else if (i == 5) {
            return notAKnotBsplineBasis.evalDx(l, i, x) + 56 * notAKnotBsplineBasis.evalDx(l, 0, x);
          } else if (i == hInv - 5) {
            return notAKnotBsplineBasis.evalDx(l, i, x) +
                   56 * notAKnotBsplineBasis.evalDx(l, hInv, x);
          } else if (i == hInv - 3) {
            return notAKnotBsplineBasis.evalDx(l, i, x) +
                   56 * notAKnotBsplineBasis.evalDx(l, hInv, x);
          } else if (i == hInv - 1) {
            return notAKnotBsplineBasis.evalDx(l, i, x) +
                   8 * notAKnotBsplineBasis.evalDx(l, hInv, x);
          } else {
            return notAKnotBsplineBasis.evalDx(l, i, x);
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
    size_t degree = getDegree();
    if ((degree != 1) && (degree != 3) && (degree != 5)) {
      throw std::runtime_error(
          "NakBsplineExtended: only B spline degrees 1, 3 and 5 are "
          "supported.");
    }

    const IT hInv = static_cast<IT>(1) << l;
    switch (getDegree()) {
        // Extension coefficients are left =[2, -1], right =[-1,2]
        // Only the ones corresponding to B-splines which are in the basis are used, i.e. indices
        // 1,end-1

      case 1:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else {
          if (i == 1) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   2 * notAKnotBsplineBasis.getIntegral(l, 0);
          } else if (i == hInv - 1) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   2 * notAKnotBsplineBasis.getIntegral(l, hInv);
          } else {
            return notAKnotBsplineBasis.getIntegral(l, i);
          }
        }

        // Constant on level 1
        // Lagrange Polynomials on level 2
        // Extended B-splines from level 3 on.

        // Extension coefficients are left =[5, -10, 10, -4], right =[-4,10,-10,5]
        // Only the ones corresponding to B-splines which are in the basis are used, i.e. indices
        // 1,3,end-3,end-1
      case 3:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else if (l == 2) {
          if (i == 1) {
            // l = 2, i = 1
            return 2.0 / 3.0;
          } else {
            // l = 2, i = 3
            return 2.0 / 3.0;
          }

        } else {
          if (i == 1) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   5 * notAKnotBsplineBasis.getIntegral(l, 0);
          } else if (i == 3) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   10 * notAKnotBsplineBasis.getIntegral(l, 0);
          } else if (i == hInv - 3) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   10 * notAKnotBsplineBasis.getIntegral(l, hInv);
          } else if (i == hInv - 1) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   5 * notAKnotBsplineBasis.getIntegral(l, hInv);
          } else {
            return notAKnotBsplineBasis.getIntegral(l, i);
          }
        }

        // Constant on level 1
        // Lagrange Polynomials on level 2
        // Extended B-splines from level 3 on.

        // Extension coefficients are left =[8,-28,42,-35,20,-6], right =[-6,20,-35,42,-28,8] on
        // level 3 and left = [8,-28,56,-70,56,-21], right = [-21,56,-70,56,-28,8] on level >= 4
        // Only the ones corresponding to B-splines which are in the basis are used, i.e. indices
        // 1,3,5,end-5,end-3,end-1
      case 5:
        if (l == 1) {
          // l = 1, i = 1
          return 1.0;
        } else if (l == 2) {
          if (i == 1) {
            // l = 2, i = 1
            return 2.0 / 3.0;
          } else {
            // l = 2, i = 3
            return 2.0 / 3.0;
          }

        } else if (l == 3) {
          if (i == 1) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   8 * notAKnotBsplineBasis.getIntegral(l, 0);
          } else if (i == 3) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   42 * notAKnotBsplineBasis.getIntegral(l, 0) +
                   20 * notAKnotBsplineBasis.getIntegral(l, hInv);
          } else if (i == 5) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   20 * notAKnotBsplineBasis.getIntegral(l, 0) +
                   42 * notAKnotBsplineBasis.getIntegral(l, hInv);
          } else if (i == 7) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   8 * notAKnotBsplineBasis.getIntegral(l, hInv);
          } else {
            return notAKnotBsplineBasis.getIntegral(l, i);
          }
        } else {
          if (i == 1) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   8 * notAKnotBsplineBasis.getIntegral(l, 0);
          } else if (i == 3) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   56 * notAKnotBsplineBasis.getIntegral(l, 0);
          } else if (i == 5) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   56 * notAKnotBsplineBasis.getIntegral(l, 0);
          } else if (i == hInv - 5) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   56 * notAKnotBsplineBasis.getIntegral(l, hInv);
          } else if (i == hInv - 3) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   56 * notAKnotBsplineBasis.getIntegral(l, hInv);
          } else if (i == hInv - 1) {
            return notAKnotBsplineBasis.getIntegral(l, i) +
                   8 * notAKnotBsplineBasis.getIntegral(l, hInv);
          } else {
            return notAKnotBsplineBasis.getIntegral(l, i);
          }
        }

      default:
        return 0.0;
    }
  }

  /**
   * Calculates the mean int b_i(x) rho(x) dx of a basis function b_i w.r.t. the probability
   * density function rho
   *
   * @param l     		level of basis function
   * @param i     		index of basis function
   * @param pdf   		probability density function
   * @param quadCoordinates coordinates of the quadrature rule to be used
   * @param quadWeights weights of the quadrature rule to be used
   * @return      		mean of basis function
   */
  inline double getMean(LT l, IT i, std::shared_ptr<sgpp::base::Distribution> pdf,
                        std::shared_ptr<sgpp::base::DataVector> quadCoordinates,
                        std::shared_ptr<sgpp::base::DataVector> quadWeights) {
    size_t degree = getDegree();
    if ((degree != 1) && (degree != 3) && (degree != 5)) {
      throw std::runtime_error(
          "NakBsplineExtended: only B spline degrees 1, 3 and 5 are "
          "supported.");
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

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override { return notAKnotBsplineBasis.getDegree(); }

 protected:
  /// B-spline basis for B-spline evaluation
  NakBsplineBasis<LT, IT> notAKnotBsplineBasis;

 private:
};

// default type-def (unsigned int for level and index)
typedef NakBsplineExtendedBasis<unsigned int, unsigned int> SNakBsplineExtendedBase;

}  // namespace base
}  // namespace sgpp
