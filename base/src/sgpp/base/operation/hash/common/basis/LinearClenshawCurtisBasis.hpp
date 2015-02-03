// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEAR_CLENSHAW_CURTIS_BASE_HPP
#define LINEAR_CLENSHAW_CURTIS_BASE_HPP

#include <cmath>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/tools/ClenshawCurtisTable.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace base {

    /**
     * Linear basis on Clenshaw-Curtis grids.
     */
    template<class LT, class IT>
    class LinearClenshawCurtisBasis: public Basis<LT, IT> {
      public:
        /**
         * @param l     level of basis function
         * @param i     index of basis function
         * @param x     evaluation point
         * @return      value of Clenshaw-Curtis linear basis function
         */
        inline double eval(LT l, IT i, double x) {
          if (l == 0) {
            // first level
            if (i == 0) {
              return 1.0 - x;
            } else {
              return x;
            }
          } else {
            // endpoints of support
            const double x0 = clenshawCurtisTable.getPoint(l, i - 1);
            const double x2 = clenshawCurtisTable.getPoint(l, i + 1);

            if ((x <= x0) || (x >= x2)) {
              // point out of support
              return 0.0;
            }

            // peak of basis function
            const double x1 = clenshawCurtisTable.getPoint(l, i);

            // linear interpolation between (x0, x1, x2), (0, 1, 0)
            if (x < x1) {
              return 1.0 - (x1 - x) / (x1 - x0);
            } else {
              return (x2 - x) / (x2 - x1);
            }
          }
        }
    };

    // default type-def (unsigned int for level and index)
    typedef LinearClenshawCurtisBasis<unsigned int, unsigned int> SLinearClenshawCurtisBase;
  }
}

#endif /* LINEAR_CLENSHAW_CURTIS_BASE_HPP */
