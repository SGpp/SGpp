/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_BASIS_BSPLINE_BOUNDARY_BSPLINEBOUNDARYBASIS_HPP
#define SGPP_BASE_BASIS_BSPLINE_BOUNDARY_BSPLINEBOUNDARYBASIS_HPP

#include <cmath>

#include "base/basis/bspline/noboundary/BsplineBasis.hpp"

namespace sg {
  namespace base {

    /**
     * B-spline basis on Boundary grids.
     */
    template <class LT, class IT>
    class BsplineBoundaryBasis {
      protected:
        /// B-spline basis for B-spline evaluation
        BsplineBasis<LT, IT> bspline_basis;

      public:
        /**
         * Default constructor.
         */
        BsplineBoundaryBasis() : bspline_basis(BsplineBasis<LT, IT>()) {
        }

        /**
         * Constructor.
         *
         * @param degree    B-spline degree, must be odd (if it's even, degree - 1 is used)
         */
        BsplineBoundaryBasis(size_t degree) : bspline_basis(BsplineBasis<LT, IT>(degree)) {
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of boundary B-spline basis function
         */
        inline double eval(LT l, IT i, double x) {
          const double hinv = static_cast<double>(1 << l);

          return bspline_basis.uniformBSpline(x * hinv - static_cast<double>(i) +
                                              static_cast<double>(bspline_basis.getDegree() + 1) / 2.0,
                                              bspline_basis.getDegree());
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of derivative of boundary B-spline basis function
         */
        inline double evalDx(LT l, IT i, double x) {
          const double hinv = static_cast<double>(1 << l);

          return hinv * bspline_basis.uniformBSplineDx(x * hinv - static_cast<double>(i) +
                 static_cast<double>(bspline_basis.getDegree() + 1) / 2.0,
                 bspline_basis.getDegree());
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of 2nd derivative of boundary B-spline basis function
         */
        inline double evalDxDx(LT l, IT i, double x) {
          const double hinv = static_cast<double>(1 << l);

          return hinv*hinv * bspline_basis.uniformBSplineDxDx(x * hinv - static_cast<double>(i) +
                 static_cast<double>(bspline_basis.getDegree() + 1) / 2.0,
                 bspline_basis.getDegree());
        }

        /**
         * @return      B-spline degree
         */
        inline size_t getDegree() const {
          return bspline_basis.degree;
        }
    };

/// typedef for standard level/index types
    typedef BsplineBoundaryBasis<unsigned int, unsigned int> SBsplineBoundaryBase;

  }
}

#endif
