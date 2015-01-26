/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALLINEARSTRETCHED_HPP
#define OPERATIONMULTIPLEEVALLINEARSTRETCHED_HPP

#include "base/operation/OperationMultipleEval.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
  namespace base {

    /**
     * This class implements OperationB for a grids with linearstretched basis ansatzfunctions without boundaries
     */
    class OperationMultipleEvalLinearStretched : public OperationMultipleEval {
      public:
        /**
         * Constructor of OperationBLinearStretched
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param dataset Pointer to the dataset that should be evaluated
         */
        OperationMultipleEvalLinearStretched(Grid &grid, DataMatrix &dataset) : OperationMultipleEval(grid, dataset) {
          this->storage = grid.getStorage();
        }

        /**
         * Destructor
         */
        virtual ~OperationMultipleEvalLinearStretched() {}

        virtual void mult(DataVector& alpha, DataVector& result);
        virtual void multTranspose(DataVector& source, DataVector& result);

      protected:
        /// Pointer to the grid's GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONMULTIPLEEVALLINEARSTRETCHED_HPP */
