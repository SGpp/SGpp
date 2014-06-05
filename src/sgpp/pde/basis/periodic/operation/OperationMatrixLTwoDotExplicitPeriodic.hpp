/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Florian Zipperle (florian.zipperle@tum.de)

#ifndef OPERATIONMATRIXLTWODOTEXPLICITPERIODIC_HPP
#define OPERATIONMATRIXLTWODOTEXPLICITPERIODIC_HPP


#include "base/operation/OperationMatrix.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/Grid.hpp"

namespace sg {
  namespace pde {

    /**
     * Explicit representation of the matrix \f$(\Phi_i,\Phi_j)_{L2}\f$ for a sparse grid
     */
    class OperationMatrixLTwoDotExplicitPeriodic: public sg::base::OperationMatrix {
      public:
        /**
         * Constructor that uses a external matrix pointer to construct the matrix,
         *
         * @param m pointer to datamatrix of size (number of grid point) x (number of grid points)
         * @param grid the sparse grid
         */
    	OperationMatrixLTwoDotExplicitPeriodic(sg::base::DataMatrix* m, sg::base::Grid* grid);
        /**
         * Constructor that creates an own matrix
         *
         * @param grid the sparse grid
         */
    	OperationMatrixLTwoDotExplicitPeriodic(sg::base::Grid* grid);

        /**
         * Destructor
         */
        virtual ~OperationMatrixLTwoDotExplicitPeriodic();

        /**
         * Implementation of standard matrix multiplication
         *
         * @param alpha DataVector that is multiplied to the matrix
         * @param result DataVector into which the result of multiplication is stored
         */
        virtual void mult(sg::base::DataVector& alpha,
                          sg::base::DataVector& result);

      private:
        /**
         * This method is used by both constructors to build the matrix
         */
        void buildMatrix(sg::base::Grid* grid);

        sg::base::DataMatrix* m_;
        bool ownsMatrix_;
    };

  } /* namespace pde */
} /* namespace sg */

#endif /* OPERATIONMATRIXLTWODOTEXPLICITPERIODIC_HPP */
