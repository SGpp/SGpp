/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
#ifndef OperationMatrixLTwoDotExplicitLinearBoundary_HPP_
#define OperationMatrixLTwoDotExplicitLinearBoundary_HPP_

#include "base/operation/OperationMatrix.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/Grid.hpp"

namespace sg {
  namespace pde {

    /**
     * Explicit representation of the matrix \f$(\Phi_i,\Phi_j)_{L2}\f$ for a sparse grid
     */
    class OperationMatrixLTwoDotExplicitLinearBoundary: public sg::base::OperationMatrix {
      public:
        /**
         * Constructor that uses a external matrix pointer to construct the matrix,
         * i.e. matrix is NOT destroyed by the destructor of OperationMatrixLTwoDotExplicitLinearBoundaryFullGrid
         *
         * @param m pointer to datamatrix of size (number of grid point) x (number of grid points)
         * @param grid the sparse grid
         */
        OperationMatrixLTwoDotExplicitLinearBoundary(sg::base::DataMatrix* m, sg::base::Grid* grid);
        /**
         * Constructor that creates an own matrix
         * i.e. matrix is destroyed by the destructor of OperationMatrixLTwoDotExplicitLinearBoundaryFullGrid
         *
         * @param grid the sparse grid
         */
        OperationMatrixLTwoDotExplicitLinearBoundary(sg::base::Grid* grid);

        /**
         * Destructor
         */
        virtual ~OperationMatrixLTwoDotExplicitLinearBoundary();

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
#endif /* OperationMatrixLTwoDotExplicitLinearBoundary_HPP_ */
