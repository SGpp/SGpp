/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_BASIS_WAVELET_MODIFIED_MODIFIEDWAVELETBASIS_HPP
#define SGPP_BASE_BASIS_WAVELET_MODIFIED_MODIFIEDWAVELETBASIS_HPP

#include <cmath>

namespace sg {
  namespace base {

    /**
     * Modified wavelet basis on Noboundary grids.
     */
    template <class LT, class IT>
    class ModWaveletBasis {
      public:
        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of modified wavelet basis function
         */
        inline double eval(LT l, IT i, double x) {
          if (l == 1) {
            // first level
            return 1.0;
          }

          const IT hinv = 1 << l;
          const double hinv_dbl = static_cast<double>(hinv);
          const double h = 1.0 / hinv_dbl;

          if ((i == 1) && (x < 1.560231504260063 * h)) {
            // left modified basis function
            return 0.501309319347014 + 1.38033323862282 * (0.560231504260063 - x * hinv_dbl + 1.0);
          } else if ((i == hinv - 1) && (x > 1.0 - 1.560231504260063 * h)) {
            // right modified basis function
            return 0.501309319347014 + 1.38033323862282 * (0.560231504260063 +
                   x * hinv_dbl - static_cast<double>(i));
          } else {
            // interior basis function
            // (or modified function, but x is in the non-modified part)
            const double t = x * hinv_dbl - static_cast<double>(i);

            if ((t > 2.0) || (t < -2.0)) {
              // out of support (cut-off)
              return 0.0;
            }

            const double t2 = t*t;
            return (1.0 - t2) * std::exp(-t2);
          }
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of derivative of modified wavelet basis function
         */
        inline double evalDx(LT l, IT i, double x) {
          if (l == 1) {
            // first level
            return 0.0;
          }

          const IT hinv = 1 << l;
          const double hinv_dbl = static_cast<double>(hinv);
          const double h = 1.0 / hinv_dbl;

          if ((i == 1) && (x < 1.560231504260063 * h)) {
            // left modified basis function
            return -1.38033323862282 * hinv_dbl;
          } else if ((i == hinv - 1) && (x > 1.0 - 1.560231504260063 * h)) {
            // right modified basis function
            return 1.38033323862282 * hinv_dbl;
          } else {
            // interior basis function
            // (or modified function, but x is in the non-modified part)
            const double t = x * hinv_dbl - static_cast<double>(i);

            if ((t > 2.0) || (t < -2.0)) {
              // out of support (cut-off)
              return 0.0;
            }

            const double t2 = t*t;
            return 2.0 * t * (t2 - 2.0) * std::exp(-t2) * hinv_dbl;
          }
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of 2nd derivative of modified wavelet basis function
         */
        inline double evalDxDx(LT l, IT i, double x) {
          if (l == 1) {
            // first level
            return 0.0;
          }

          const IT hinv = 1 << l;
          const double hinv_dbl = static_cast<double>(hinv);
          const double h = 1.0 / hinv_dbl;

          if ((i == 1) && (x < 1.560231504260063 * h)) {
            // left modified basis function
            return 0.0;
          } else if ((i == hinv - 1) && (x > 1.0 - 1.560231504260063 * h)) {
            // right modified basis function
            return 0.0;
          } else {
            // interior basis function
            // (or modified function, but x is in the non-modified part)
            const double t = x * hinv_dbl - static_cast<double>(i);

            if ((t > 2.0) || (t < -2.0)) {
              // out of support (cut-off)
              return 0.0;
            }

            const double t2 = t*t;
            const double t4 = t2*t2;
            return -2.0 * (2.0*t4 - 7.0*t2 + 2.0) * std::exp(-t2) * hinv_dbl*hinv_dbl;
          }
        }
    };

/// typedef for standard level/index types
    typedef ModWaveletBasis<unsigned int, unsigned int> SModWaveletBase;

  }
}

#endif
