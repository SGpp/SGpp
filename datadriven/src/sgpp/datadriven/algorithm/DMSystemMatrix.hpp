/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DMSYSTEMMATRIX_HPP
#define DMSYSTEMMATRIX_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/Grid.hpp"
#include "base/operation/OperationMultipleEval.hpp"
#include "base/operation/OperationMatrix.hpp"

#include "datadriven/algorithm/DMSystemMatrixBase.hpp"

namespace sg {
  namespace datadriven {

    /**
     * Class that implements the virtual class sg::base::OperationMatrix for the
     * application of classification for the Systemmatrix
     */
    class DMSystemMatrix : public DMSystemMatrixBase {
      private:
        /// sg::base::OperationMatrix, the regularisation mehtod
        sg::base::OperationMatrix* C;
        /// OperationB for calculating the data matrix
        //OperationMultiEval* B;
        sg::base::OperationMultipleEval* B;

        sg::base::Grid &grid;

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param trainData reference to sg::base::DataVector that contains the training data
         * @param C the regression functional
         * @param lambda the lambda, the regression parameter
         */
        DMSystemMatrix(sg::base::Grid& grid, sg::base::DataMatrix& trainData, sg::base::OperationMatrix& C, double lambda);

        /**
         * Std-Destructor
         */
        virtual ~DMSystemMatrix();

        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

        /**
         * Generates the right hand side of the classification equation
         *
         * @param classes the class information of the training data
         * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
         */
        virtual void generateb(sg::base::DataVector& classes, sg::base::DataVector& b);
    };

  }
}

#endif /* DMSYSTEMMATRIX_HPP */
