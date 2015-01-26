/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Florian Zipperle (florian.zipperle@tum.de)

#ifndef OPERATIONMULTIPLEEVALPERIODIC_HPP
#define OPERATIONMULTIPLEEVALPERIODIC_HPP

#include <sgpp/base/operation/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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
    	OperationMultipleEvalPeriodic(Grid &grid, DataMatrix &dataset) : OperationMultipleEval(grid, dataset) {
          this->storage = grid.getStorage();
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
