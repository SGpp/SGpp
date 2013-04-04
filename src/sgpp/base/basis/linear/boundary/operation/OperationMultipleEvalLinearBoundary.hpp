/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALLINEARBOUNDARY_HPP
#define OPERATIONMULTIPLEEVALLINEARBOUNDARY_HPP

#include "base/operation/OperationMultipleEval.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
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
        OperationMultipleEvalLinearBoundary(GridStorage* storage, DataMatrix* dataset) : OperationMultipleEval(dataset) {
          this->storage = storage;
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
