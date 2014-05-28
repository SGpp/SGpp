/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Florian Zipperle (florian.zipperle@tum.de)

#ifndef OPERATIONMULTIPLEEVALPERIODIC_HPP
#define OPERATIONMULTIPLEEVALPERIODIC_HPP

#include "base/operation/OperationMultipleEval.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
  namespace base {

    /**
     * This class implements OperationMultipleEval for a grids with periodic linear basis ansatzfunctions
     *
     * @version $HEAD$
     */
    class OperationMultipleEvalPeriodic : public OperationMultipleEval {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GridStorage object
         * @param dataset the dataset that should be evaluated
         */
    	OperationMultipleEvalPeriodic(GridStorage* storage, DataMatrix* dataset) : OperationMultipleEval(dataset) {
          this->storage = storage;
        }

        /**
         * Destructor
         */
        virtual ~OperationMultipleEvalPeriodic() {}

        virtual void mult(DataVector& alpha, DataVector& result);
        virtual void multTranspose(DataVector& source, DataVector& result);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONMULTIPLEEVALPERIODIC_HPP */
