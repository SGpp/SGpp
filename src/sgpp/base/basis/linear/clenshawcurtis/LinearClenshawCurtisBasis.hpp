/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_BASIS_LINEAR_CLENSHAWCURTIS_LINEARCLENSHAWCURTISBASIS_HPP
#define SGPP_BASE_BASIS_LINEAR_CLENSHAWCURTIS_LINEARCLENSHAWCURTISBASIS_HPP

#include "base/tools/CosineTable.hpp"

#include <cmath>

namespace sg {
  namespace base {

    /**
     * Linear basis on Clenshaw-Curtis grids.
     */
    template <class LT, class IT>
    class LinearClenshawCurtisBasis {
      protected:
        const CosineTable* cosine_table;

      public:
        /**
         * Constructor.
         *
         * @param cosine_table  cosine table for faster cosine evaluation (optional)
         */
        LinearClenshawCurtisBasis(const CosineTable* cosine_table = NULL) :
          cosine_table(cosine_table) {
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
            const double h = 1.0 / static_cast<double>(1 << l);
            // endpoints of support
            const double x0 = clenshawCurtisPoint(h, i - 1);
            const double x2 = clenshawCurtisPoint(h, i + 1);

            if ((x <= x0) || (x >= x2)) {
              // point out of support
              return 0.0;
            }

            // peak of basis function
            const double x1 = clenshawCurtisPoint(h, i);

            // linear interpolation between (x0, x1, x2), (0, 1, 0)
            if (x < x1) {
              return 1.0 - (x1 - x) / (x1 - x0);
            } else {
              return (x2 - x) / (x2 - x1);
            }
          }
        }
    };

/// typedef for standard level/index types
    typedef LinearClenshawCurtisBasis<unsigned int, unsigned int> SLinearClenshawCurtisBase;

  }
}

#endif
