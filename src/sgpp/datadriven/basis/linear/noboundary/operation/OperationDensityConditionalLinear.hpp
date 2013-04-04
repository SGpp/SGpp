/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Benjamin Peherstorfer (pehersto@in.tum.de)

#ifndef OPERATIONDENSITYCONDITIONALLINEAR_HPP
#define OPERATIONDENSITYCONDITIONALLINEAR_HPP

#include "base/grid/Grid.hpp"
#include "datadriven/operation/OperationDensityConditional.hpp"
#include <cstring>

namespace sg {
  namespace datadriven {

    /**
     * Marginalize Probability Density Function
     */

    class OperationDensityConditionalLinear : public OperationDensityConditional {
      public:
        OperationDensityConditionalLinear(base::Grid* grid) : grid(grid) {}
        virtual ~OperationDensityConditionalLinear() {}

        /**
         * Marginalizes (Density) Functions
         *
         * @param alpha Coefficient vector for current grid
         * @param mg Referenz of grid pointer
         * @param malpha Coefficient vector for new grid (mg). Will be resized.
         * @param mdim Marginalize in dimension mdim
         * @param xbar Point at which to conditionalize
         */
        void doConditional(base::DataVector& alpha, base::Grid*& mg, base::DataVector& malpha, unsigned int mdim, double xbar);

      protected:
        base::Grid* grid;
    };

  }
}
#endif /* OPERATIONDENSITYCONDITIONALLINEAR_HPP */
