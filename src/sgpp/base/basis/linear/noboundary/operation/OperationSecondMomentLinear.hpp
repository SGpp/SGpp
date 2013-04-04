/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de)
// @author Benjamin

#ifndef OPERATIONSECONDMOMENTLINEAR_HPP
#define OPERATIONSECONDMOMENTLINEAR_HPP

#include "base/operation/OperationSecondMoment.hpp"
#include "base/grid/Grid.hpp"

namespace sg {
  namespace base {

    /**
     * SecondMomemnt of sparse grid function, linear grid without boundaries
     */
    class OperationSecondMomentLinear : public OperationSecondMoment {
      public:
        /**
         * Constructor of OperationSecondMomentLinear
         *
         * @param storage Pointer to the grid's GridStorage object
         */
        OperationSecondMomentLinear(GridStorage* storage) : storage(storage) {}

        virtual ~OperationSecondMomentLinear() {}

        /**
         * Compute second moment of the function
         * @f[ \int_{\Omega} x^2\cdot f(x) dx. @f]
         *
         * @param alpha Coefficient vector for current grid
         */
        virtual double doQuadrature(DataVector& alpha);

      protected:
        // Pointer to the grid's GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONSECONDMOMENTLINEAR_HPP */
