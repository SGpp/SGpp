/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_BASIS_BSPLINE_CLENSHAWCURTIS_BSPLINECLENSHAWCURTISBASIS_HPP
#define SGPP_BASE_BASIS_BSPLINE_CLENSHAWCURTIS_BSPLINECLENSHAWCURTISBASIS_HPP

#include <cmath>

#include "base/basis/bspline/noboundary/BsplineBasis.hpp"
#include "base/tools/CosineTable.hpp"

namespace sg {
  namespace base {

    /**
     * B-spline basis on Clenshaw-Curtis grids.
     */
    template <class LT, class IT>
    class BsplineClenshawCurtisBasis {
      protected:
        /// B-spline basis for B-spline evaluation
        BsplineBasis<LT, IT> bspline_basis;
        const CosineTable* cosine_table;
        std::vector<double> knots;

      public:
        /**
         * Default constructor.
         */
        BsplineClenshawCurtisBasis() : bspline_basis(BsplineBasis<LT, IT>()), cosine_table(NULL) {
        }

        /**
         * Constructor.
         *
         * @param degree        B-spline degree, must be odd (if it's even, degree - 1 is used)
         * @param cosine_table  cosine table for faster cosine evaluation (optional)
         */
        BsplineClenshawCurtisBasis(size_t degree, const CosineTable* cosine_table = NULL) :
          bspline_basis(BsplineBasis<LT, IT>(degree)),
          cosine_table(cosine_table),
          knots(std::vector<double>(degree+2, 0.0)) {
        }

        /**
         * @param x     evaluation point
         * @param p     B-spline degree
         * @param k     index of B-spline in the knot sequence
         * @param xi    knot sequence
         * @return      value of non-uniform B-spline with knots
         *              \f$\{\xi_k, ... \xi_{k+p+1}\}\f$
         */
        inline double nonUniformBSpline(double x, size_t p, size_t k,
                                        const std::vector<double>& xi) const {
          if (p == 0) {
            // characteristic function of [xi[k], xi[k+1])
            return (((x >= xi[k]) && (x < xi[k+1])) ? 1.0 : 0.0);
          } else if ((x < xi[k]) || (x >= xi[k+p+1])) {
            // out of support
            return 0.0;
          } else {
            // Cox-de-Boor recursion
            return (x - xi[k]) / (xi[k+p] - xi[k]) * nonUniformBSpline(x, p-1, k, xi)
                   + (1.0 - (x - xi[k+1]) / (xi[k+p+1] - xi[k+1]))
                   * nonUniformBSpline(x, p-1, k+1, xi);
          }
        }

        /**
         * @param x     evaluation point
         * @param p     B-spline degree
         * @param k     index of B-spline in the knot sequence
         * @param xi    knot sequence
         * @return      value of derivative of non-uniform B-spline with knots
         *              \f$\{\xi_k, ... \xi_{k+p+1}\}\f$
         */
        inline double nonUniformBSplineDx(double x, size_t p, size_t k,
                                          const std::vector<double>& xi) const {
          if (p == 0) {
            return 0.0;
          } else if ((x < xi[k]) || (x >= xi[k+p+1])) {
            return 0.0;
          } else {
            const double p_double = static_cast<double>(p);

            return p_double / (xi[k+p] - xi[k]) * nonUniformBSpline(x, p-1, k, xi)
                   - p_double / (xi[k+p+1] - xi[k+1]) * nonUniformBSpline(x, p-1, k+1, xi);
          }
        }

        /**
         * @param x     evaluation point
         * @param p     B-spline degree
         * @param k     index of B-spline in the knot sequence
         * @param xi    knot sequence
         * @return      value of 2nd derivative of non-uniform B-spline with knots
         *              \f$\{\xi_k, ... \xi_{k+p+1}\}\f$
         */
        inline double nonUniformBSplineDxDx(double x, size_t p, size_t k,
                                            const std::vector<double>& xi) const {
          if (p <= 1) {
            return 0.0;
          } else if ((x < xi[k]) || (x >= xi[k+p+1])) {
            return 0.0;
          } else {
            const double p_double = static_cast<double>(p);
            const double alpha_k_p = p_double / (xi[k+p] - xi[k]);
            const double alpha_kp1_p = p_double / (xi[k+p+1] - xi[k+1]);
            const double alpha_k_pm1 = (p_double - 1.0) / (xi[k+p-1] - xi[k]);
            const double alpha_kp1_pm1 = (p_double - 1.0) / (xi[k+p] - xi[k+1]);
            const double alpha_kp2_pm1 = (p_double - 1.0) / (xi[k+p+1] - xi[k+2]);

            return alpha_k_p * alpha_k_pm1 * nonUniformBSpline(x, p-2, k, xi)
                   - (alpha_k_p + alpha_kp1_p) * alpha_kp1_pm1 *
                   nonUniformBSpline(x, p-2, k+1, xi)
                   + alpha_kp1_p * alpha_kp2_pm1 * nonUniformBSpline(x, p-2, k+2, xi);
          }
        }

        /**
         * @param h     grid width (= 2^(-l)) of the grid point
         * @param i     index of the grid point
         * @return      i-th Clenshaw-Curtis grid point with level -log2(h)
         */
        inline double clenshawCurtisPoint(double h, IT i) const {
          if (cosine_table == NULL) {
            return (cos(M_PI * (1.0 - static_cast<double>(i) * h)) + 1.0) / 2.0;
          } else {
            return (cosine_table->lookUp(
                      M_PI * (1.0 - static_cast<double>(i) * h)) + 1.0) / 2.0;
          }
        }

        /**
         * Construct the (p+2) Clenshaw-Curtis knots of a B-spline basis function
         * and save them in knots.
         *
         * @param l     level of basis function
         * @param i     index of basis function
         */
        inline void constructKnots(LT l, IT i) {
          const IT hinv = 1 << l;
          const double h = 1.0 / static_cast<double>(hinv);
          const size_t p = bspline_basis.getDegree();

          knots[(p+1)/2] = clenshawCurtisPoint(h, i);

          if (i < (p+1)/2) {
            // grid point index is too far on the left
            // ==> extrapolate grid points linearly
            size_t a = (p+1)/2 - i;

            for (size_t j = a; j < (p+1)/2; j++) {
              knots[j] = clenshawCurtisPoint(h, static_cast<IT>(j - a));
            }

            double h = knots[a+1] - knots[a];

            // equivalent to "for (int j = a-1; j >= 0; j--)"
            for (size_t j = a; j-- > 0; ) {
              knots[j] = knots[j+1] - h;
            }
          } else {
            // all grid points on the left can be calculated
            for (size_t j = 0; j < (p+1)/2; j++) {
              knots[j] = clenshawCurtisPoint(h, static_cast<IT>(i - (p+1)/2 + j));
            }
          }

          if (i + (p+1)/2 > hinv) {
            // grid point index is too far on the right
            // ==> extrapolate grid points linearly
            size_t b = hinv + (p+1)/2 - i;

            for (size_t j = (p+1)/2 + 1; j <= b; j++) {
              knots[j] = clenshawCurtisPoint(h, static_cast<IT>(i - (p+1)/2 + j));
            }

            double h = knots[b] - knots[b-1];

            for (size_t j = b+1; j < p+2; j++) {
              knots[j] = knots[j-1] + h;
            }
          } else {
            // all grid points on the right can be calculated
            for (size_t j = (p+1)/2 + 1; j < p+2; j++) {
              knots[j] = clenshawCurtisPoint(h, static_cast<IT>(i - (p+1)/2 + j));
            }
          }
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of Clenshaw-Curtis B-spline basis function
         */
        inline double eval(LT l, IT i, double x) {
          if (l == 0) {
            return bspline_basis.uniformBSpline(x - static_cast<double>(i) +
                                                static_cast<double>(bspline_basis.getDegree() + 1) / 2.0,
                                                bspline_basis.getDegree());
          } else {
            constructKnots(l, i);
            return nonUniformBSpline(x, bspline_basis.getDegree(), 0, knots);
          }
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of derivative of Clenshaw-Curtis B-spline basis function
         */
        inline double evalDx(LT l, IT i, double x) {
          if (l == 0) {
            return bspline_basis.uniformBSplineDx(x - static_cast<double>(i) +
                                                  static_cast<double>(bspline_basis.getDegree() + 1) / 2.0,
                                                  bspline_basis.getDegree());
          } else {
            constructKnots(l, i);
            return nonUniformBSplineDx(x, bspline_basis.getDegree(), 0, knots);
          }
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of 2nd derivative of Clenshaw-Curtis B-spline basis function
         */
        inline double evalDxDx(LT l, IT i, double x) {
          if (l == 0) {
            return bspline_basis.uniformBSplineDxDx(x - static_cast<double>(i) +
                                                    static_cast<double>(bspline_basis.getDegree() + 1) / 2.0,
                                                    bspline_basis.getDegree());
          } else {
            constructKnots(l, i);
            return nonUniformBSplineDxDx(x, bspline_basis.getDegree(), 0, knots);
          }
        }

        /**
         * @return      B-spline degree
         */
        inline size_t getDegree() const {
          return bspline_basis.degree;
        }
    };

/// typedef for standard level/index types
    typedef BsplineClenshawCurtisBasis<unsigned int, unsigned int> SBsplineClenshawCurtisBase;

  }
}

#endif
