// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARSTRETCHEDBOUNDARYBASE_HPP
#define LINEARSTRETCHEDBOUNDARYBASE_HPP

#include <cmath>
#include <sgpp/base/grid/common/Stretching.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * linearstretched basis functions with boundaries
     * And here we have another implicit dependence on tensor products
     *
     */
    template<class LT, class IT>
    class LinearStretchedBoundaryBasis: public LinearBoundaryBasis<LT, IT> {
      public:
        /*  *
         * Evaluate a basis function.
         * Has a dependence on the absolute position of grid point and support.
         *
         * @param level the level of the current basis function
         * @param index the index of the current basis function
         * @param p the absolute position of the evaluation point

        float_t eval(LT level, IT index, float_t p)
        {
          if (level == 0)
          {
            if (index == 0)
            {
              return 1.0 - p;
            }
            if (index == 1)
            {
              return p;
            }
          }
          else
          {
            return 1.0 - fabs((1<<level) * p - index);
          }
          // should not happen
          return 0.0;
        }

         *
         * Evaluate a basis function with an offset and scaling factor
         * Has a dependence on the absolute position of grid point and support.
         *
         * @param level the level of the current basis function
         * @param index the index of the current basis function
         * @param p the absolute position of the evaluation point
         * @param q the scaling factor of the basis function
         * @param t the offset of the basis function

        float_t eval(LT level, IT index, float_t p, float_t q, float_t t)
        {
          if (level == 0)
          {
            if (index == 0)
            {
              return 1.0 - ((1.0/q)*(p-t));
            }
            if (index == 1)
            {
              return ((1.0/q)*(p-t));
            }
          }
          else
          {
            return 1.0 - ((1.0/q)*(fabs(((1<<level)*(p-t))-(q*index))));
          }
          // should not happen
          return 0.0;
        }*/

        /*
         * evaluate a basis function
         * Has a dependence on the position of two grid points with values 1 and 0 and the
         * support position
         */
        float_t eval(float_t p, float_t pos0, float_t pos1) {
          return (p - pos0) / (pos1 - pos0);
        }

        ///TODO: Index and level is not necessary, maybe function could be changed
        /*
         * evaluate a basis function
         * Has a dependence on the position of two grid points with values 1 and 0 and the
         * support position
         */
        float_t eval(LT level, IT index, float_t p, float_t pos0, float_t pos1) {
          //    if(level == 0){
          //      if(index == 0){
          //
          //      }
          //      if(index == 1){
          //
          //      }
          //
          //    }
          //    else{
          //      return (p-pos0)/(pos1-pos0);
          //    }
          //    ///Should not happen
          //    return 0.0;
          return (p - pos0) / (pos1 - pos0);
        }

    };

    // default type-def (unsigned int for level and index)
    typedef LinearStretchedBoundaryBasis<unsigned int, unsigned int> SLinearStretchedBoundaryBase;
  }
}

#endif /* LINEARSTRETCHEDBOUNDARYBASE_HPP */