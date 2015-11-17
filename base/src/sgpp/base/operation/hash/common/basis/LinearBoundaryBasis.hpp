// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEAR_BOUNDARY_BASE_HPP
#define LINEAR_BOUNDARY_BASE_HPP

#include <algorithm>
#include <cmath>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace base {

    /**
     * Linear basis on Boundary grids.
     */
    template<class LT, class IT>
    class LinearBoundaryBasis: public Basis<LT, IT> {
      public:
        /**
         * Destructor.
         */
        virtual ~LinearBoundaryBasis() override {
        }

        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of boundary linear basis function
         */
        inline virtual float_t eval(LT l, IT i, float_t x) override {
          if (l == 0) {
            // first level
            if (i == 0) {
              return 1.0 - x;
            } else {
              return x;
            }
          } else {
            return std::max(1.0 - std::abs(static_cast<float_t>(static_cast<IT>(1) << l) * x -
                                           static_cast<float_t>(i)), 0.0);
          }
        }

        /**
         * Evaluate basis function with offset and scaling factor.
         *
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @param q     scaling factor of basis function
         * @param t     offset of basis function
         */
        inline virtual float_t eval(LT l, IT i, float_t x, float_t q, float_t t) {
          return eval(l, i, (x - t) / q);
        }
    };

    // default type-def (unsigned int for level and index)
    typedef LinearBoundaryBasis<unsigned int, unsigned int> SLinearBoundaryBase;

  }
}

#endif /* LINEAR_BOUNDARY_BASE_HPP */
