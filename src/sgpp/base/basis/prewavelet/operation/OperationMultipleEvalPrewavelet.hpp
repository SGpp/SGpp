/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de)
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALPREWAVELET_HPP
#define OPERATIONMULTIPLEEVALPREWAVELET_HPP

#include "base/operation/OperationMultipleEval.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
  namespace base {

    /**
     * This class implements OperationMultipleEval for a grids with prewavelet ansatzfunctions without boundaries
     */
    class OperationMultipleEvalPrewavelet : public OperationMultipleEval {
      public:
        /**
         * Constructor of OperationMultipleEvalPrewavelet
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param dataset Dataset
         */
        OperationMultipleEvalPrewavelet(Grid &grid, DataMatrix &dataset) : OperationMultipleEval(grid, dataset) {
          this->storage = grid.getStorage();
        }

        /**
         * Destructor
         */
        virtual ~OperationMultipleEvalPrewavelet() {}

        virtual void mult(DataVector& alpha, DataVector& result);
        virtual void multTranspose(DataVector& source, DataVector& result);

      protected:
        /// Pointer to the grid's GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONMULTIPLEEVALPREWAVELET_HPP */
