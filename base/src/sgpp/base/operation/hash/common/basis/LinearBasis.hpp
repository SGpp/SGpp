// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEAR_BASE_HPP
#define LINEAR_BASE_HPP

#include <cmath>
#include <algorithm>

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace base {

    /**
     * Linear basis on Noboundary grids.
     */
    template<class LT, class IT>
    class LinearBasis: public Basis<LT, IT> {
      public:
        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of linear basis function
         */
        inline double eval(LT l, IT i, double x) {
          return std::max(1.0 - std::abs(static_cast<double>(static_cast<IT>(1) << l) * x -
                                         static_cast<double>(i)), 0.0);
        }
    };

    // default type-def (unsigned int for level and index)
    typedef LinearBasis<unsigned int, unsigned int> SLinearBase;

  }
}

#endif /* LINEAR_BASE_HPP */
