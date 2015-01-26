// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALPOLY_HPP
#define OPERATIONMULTIPLEEVALPOLY_HPP

#include <sgpp/base/operation/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class implements OperationMultipleEval for a grids with poly basis ansatzfunctions
     */
    class OperationMultipleEvalPoly : public OperationMultipleEval {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GridStorage object
         * @param degree the polynom's max. degree
         * @param dataset Dataset
         */
        OperationMultipleEvalPoly(Grid &grid, size_t degree, DataMatrix &dataset) : OperationMultipleEval(grid, dataset), base(degree) {
          this->storage = storage;
        }

        /**
         * Destructor
         */
        virtual ~OperationMultipleEvalPoly() {}

        virtual void mult(DataVector& alpha, DataVector& result);
        virtual void multTranspose(DataVector& source, DataVector& result);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
        /// Poly Basis object
        SPolyBase base;
    };

  }
}

#endif /* OPERATIONMULTIPLEEVALPOLY_HPP */