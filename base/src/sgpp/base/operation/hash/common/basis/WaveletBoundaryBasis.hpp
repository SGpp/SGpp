// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WAVELET_BOUNDARY_BASE_HPP
#define WAVELET_BOUNDARY_BASE_HPP

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
         * Destructor.
         */
        virtual ~WaveletBoundaryBasis() override {
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of boundary wavelet basis function
         */
        inline virtual float_t eval(LT l, IT i, float_t x) override {
          const float_t hinv = static_cast<float_t>(1 << l);
          const float_t t = x * hinv - static_cast<float_t>(i);

          if ((t >= 2.0) || (t <= -2.0)) {
            // out of support (cut-off)
            return 0.0;
          }

          const float_t t2 = t * t;
          return (1.0 - t2) * std::exp(-t2);
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of derivative of boundary wavelet basis function
         */
        inline float_t evalDx(LT l, IT i, float_t x) {
          const float_t hinv = static_cast<float_t>(1 << l);
          const float_t t = x * hinv - static_cast<float_t>(i);

          if ((t >= 2.0) || (t <= -2.0)) {
            // out of support (cut-off)
            return 0.0;
          }

          const float_t t2 = t * t;
          return 2.0 * t * (t2 - 2.0) * std::exp(-t2) * hinv;
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of 2nd derivative of boundary wavelet basis function
         */
        inline float_t evalDxDx(LT l, IT i, float_t x) {
          const float_t hinv = static_cast<float_t>(1 << l);
          const float_t t = x * hinv - static_cast<float_t>(i);

          if ((t >= 2.0) || (t <= -2.0)) {
            // out of support (cut-off)
            return 0.0;
          }

          const float_t t2 = t * t;
          return -2.0 * (2.0 * t2 * t2 - 7.0 * t2 + 2.0) * std::exp(-t2) * hinv * hinv;
        }
    };

    // default type-def (unsigned int for level and index)
    typedef WaveletBoundaryBasis<unsigned int, unsigned int> SWaveletBoundaryBase;
  }
}

#endif /* WAVELET_BOUNDARY_BASE_HPP */
