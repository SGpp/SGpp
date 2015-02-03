// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WAVELET_BASE_HPP
#define WAVELET_BASE_HPP

#include <cmath>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace base {

    /**
     * Wavelet basis on Boundary grids.
     */
    template <class LT, class IT>
    class WaveletBoundaryBasis: public Basis<LT, IT> {
      public:
        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of boundary wavelet basis function
         */
        inline double eval(LT l, IT i, double x) {
          const double hinv = static_cast<double>(1 << l);
          const double t = x * hinv - static_cast<double>(i);

          if ((t >= 2.0) || (t <= -2.0)) {
            // out of support (cut-off)
            return 0.0;
          }

          const double t2 = t * t;
          return (1.0 - t2) * std::exp(-t2);
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of derivative of boundary wavelet basis function
         */
        inline double evalDx(LT l, IT i, double x) {
          const double hinv = static_cast<double>(1 << l);
          const double t = x * hinv - static_cast<double>(i);

          if ((t >= 2.0) || (t <= -2.0)) {
            // out of support (cut-off)
            return 0.0;
          }

          const double t2 = t * t;
          return 2.0 * t * (t2 - 2.0) * std::exp(-t2) * hinv;
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of 2nd derivative of boundary wavelet basis function
         */
        inline double evalDxDx(LT l, IT i, double x) {
          const double hinv = static_cast<double>(1 << l);
          const double t = x * hinv - static_cast<double>(i);

          if ((t >= 2.0) || (t <= -2.0)) {
            // out of support (cut-off)
            return 0.0;
          }

          const double t2 = t * t;
          return -2.0 * (2.0 * t2 * t2 - 7.0 * t2 + 2.0) * std::exp(-t2) * hinv * hinv;
        }
    };

    // default type-def (unsigned int for level and index)
    typedef WaveletBoundaryBasis<unsigned int, unsigned int> SWaveletBoundaryBase;
  }
}

#endif /* WAVELET_BASE_HPP */
