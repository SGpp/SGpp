// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONFUNDAMENTALSPLINE_HPP
#define OPERATIONHIERARCHISATIONFUNDAMENTALSPLINE_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/type/FundamentalSplineGrid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Hierarchisation on sparse grid, fundamental spline basis
     */
    class OperationHierarchisationFundamentalSpline :
      public OperationHierarchisation {
      public:
        /**
         * Constructor of OperationHierarchisationFundamentalSpline
         *
         * @param grid Pointer to the grid
         */
        OperationHierarchisationFundamentalSpline(FundamentalSplineGrid* grid);

        /**
         * Destructor.
         */
        virtual ~OperationHierarchisationFundamentalSpline();

        virtual void doHierarchisation(DataVector& node_values);
        virtual void doDehierarchisation(DataVector& alpha);

      protected:
        /// grid
        FundamentalSplineGrid* grid;
    };

  }
}

#endif /* OPERATIONHIERARCHISATIONFUNDAMENTALSPLINE_HPP */
