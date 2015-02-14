// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONDENSITYCONDITIONALLINEAR_HPP
#define OPERATIONDENSITYCONDITIONALLINEAR_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityConditional.hpp>
#include <cstring>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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
        void doConditional(base::DataVector& alpha, base::Grid*& mg, base::DataVector& malpha, unsigned int mdim, float_t xbar);

      protected:
        base::Grid* grid;
    };

  }
}
#endif /* OPERATIONDENSITYCONDITIONALLINEAR_HPP */