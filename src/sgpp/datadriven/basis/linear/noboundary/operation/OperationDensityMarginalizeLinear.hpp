/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Benjamin Peherstorfer (pehersto@in.tum.de)

#ifndef OPERATIONDENSITYMARGINALIZELINEAR_HPP
#define OPERATIONDENSITYMARGINALIZELINEAR_HPP

#include "base/grid/Grid.hpp"
#include "datadriven/operation/OperationDensityMarginalize.hpp"
#include <cstring>

namespace sg {
  namespace datadriven {

    /**
     * Marginalize Probability Density Function
     */

    class OperationDensityMarginalizeLinear : public OperationDensityMarginalize {
      public:
        OperationDensityMarginalizeLinear(base::Grid* grid) : grid(grid) {}
        virtual ~OperationDensityMarginalizeLinear() {}

        /**
         * Marginalizes (Density) Functions
         *
         * @param alpha Coefficient vector for current grid
         * @param mg Referenz of grid pointer
         * @param malpha Coefficient vector for new grid (mg). Will be resized.
         * @param mdim Marginalize in dimension mdim
         */
        void doMarginalize(base::DataVector& alpha, base::Grid*& mg, base::DataVector& malpha, unsigned int mdim);

      protected:
        base::Grid* grid;
    };

  }
}
#endif /* OPERATIONDENSITYMARGINALIZELINEAR_HPP */
