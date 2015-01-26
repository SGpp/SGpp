/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALLINEARBOUNDARY_HPP
#define OPERATIONMULTIPLEEVALLINEARBOUNDARY_HPP

#include <sgpp/base/operation/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class implements OperationMultipleEval for a grids with linear basis ansatzfunctions
     * with boundaries
     *
     * @version $HEAD$
     */
    class OperationMultipleEvalLinearBoundary : public OperationMultipleEval {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GirdStorage object
         * @param dataset the dataset the should be evaluated
         */
        OperationMultipleEvalLinearBoundary(Grid &grid, DataMatrix &dataset) : OperationMultipleEval(grid, dataset) {
          this->storage = grid.getStorage();
        }

        /**
         * Destructor
         */
        virtual ~OperationMultipleEvalLinearBoundary() {}

        virtual void mult(DataVector& alpha, DataVector& result);
        virtual void multTranspose(DataVector& source, DataVector& result);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONMULTIPLEEVALLINEARBOUNDARY_HPP */
