// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARPERIODICBASE_HPP
#define LINEARPERIODICBASE_HPP

#include <cmath>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * linear basis functions with boundaries
     * And here we have another implicit dependence on tensor products
     *
     */
    template<class LT, class IT>
    class LinearPeriodicBasis: public Basis<LT, IT>  {
      public:
        /**
         * Evaluate a basis function.
         * Has a dependence on the absolute position of grid point and support.
         *
         * @param level the level of the current basis function
         * @param index the index of the current basis function
         * @param p the absolute position of the evaluation point
         */
        float_t eval(LT level, IT index, float_t p) {
          if (level == 0) {
            return fabs(2 * p - 1);
          } else {
            return 1.0 - fabs((1 << level) * p - index);
          }

          // should not happen
          return 0.0;
        }

        /**
         * Evaluate a basis function with an offset and scaling factor
         * Has a dependence on the absolute position of grid point and support.
         *
         * @param level the level of the current basis function
         * @param index the index of the current basis function
         * @param p the absolute position of the evaluation point
         * @param q the scaling factor of the basis function
         * @param t the offset of the basis function
         */
        float_t eval(LT level, IT index, float_t p, float_t q, float_t t) {
          if (level == 0) {
            return ((1.0 / q) * (fabs((2 * (p - t)) - (q))));
          } else {
            return 1.0 - ((1.0 / q) * (fabs(((1 << level) * (p - t)) - (q * index))));
          }

          // should not happen
          return 0.0;
        }
    };

    // default type-def (unsigned int for level and index)
    typedef LinearPeriodicBasis<unsigned int, unsigned int> SLinearPeriodicBasis;

  }
}

#endif /* LINEARPERIODICBASE_HPP */