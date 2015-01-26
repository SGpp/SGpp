/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Florian Zipperle (florian.zipperle@tum.de)

#ifndef OPERATIONMATRIXLTWODOTPERIODIC_HPP
#define OPERATIONMATRIXLTWODOTPERIODIC_HPP

#include <sgpp/base/operation/OperationMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>

namespace sg {
  namespace pde {

    /**
     * Implements the standard L 2 scalar product on periodic grids
     *
     * @version $HEAD$
     */
    class OperationMatrixLTwoDotPeriodic: public sg::base::OperationMatrix {
      public:
        /**
         * Constructor
         *
         * @param gird a referenz to the grid
         */
    	OperationMatrixLTwoDotPeriodic(sg::base::GridStorage* gridStorage);

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
		virtual void mult(sg::base::DataVector& alpha,
						  sg::base::DataVector& result);
      protected:

		sg::base::GridStorage* gridStorage;
    };

  }
}


#endif /* OPERATIONMATRIXLTWODOTPERIODIC_HPP */
