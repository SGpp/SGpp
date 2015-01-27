// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONMATRIXLTWODOTPERIODIC_HPP
#define OPERATIONMATRIXLTWODOTPERIODIC_HPP

#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    /**
     * Implements the standard L 2 scalar product on periodic grids
     *
     * @version $HEAD$
     */
    class OperationMatrixLTwoDotPeriodic: public SGPP::base::OperationMatrix {
      public:
        /**
         * Constructor
         *
         * @param gird a referenz to the grid
         */
    	OperationMatrixLTwoDotPeriodic(SGPP::base::GridStorage* gridStorage);

        /**
         * Destructor
         */
        virtual ~OperationMatrixLTwoDotPeriodic();

        /**
		 * Implementation of standard matrix multiplication
		 *
		 * @param alpha DataVector that is multiplied to the matrix
		 * @param result DataVector into which the result of multiplication is stored
		 */
		virtual void mult(SGPP::base::DataVector& alpha,
						  SGPP::base::DataVector& result);
      protected:

		SGPP::base::GridStorage* gridStorage;
    };

  }
}


#endif /* OPERATIONMATRIXLTWODOTPERIODIC_HPP */