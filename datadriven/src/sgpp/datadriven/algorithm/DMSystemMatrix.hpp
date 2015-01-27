// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIX_HPP
#define DMSYSTEMMATRIX_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * Class that implements the virtual class SGPP::base::OperationMatrix for the
     * application of classification for the Systemmatrix
     */
    class DMSystemMatrix : public DMSystemMatrixBase {
      private:
        /// SGPP::base::OperationMatrix, the regularisation mehtod
        SGPP::base::OperationMatrix* C;
        /// OperationB for calculating the data matrix
        //OperationMultiEval* B;
        SGPP::base::OperationMultipleEval* B;

        SGPP::base::Grid& grid;

      public:
        /**
         * Std-Constructor
         *
         * @param grid reference to the sparse grid
         * @param trainData reference to SGPP::base::DataVector that contains the training data
         * @param C the regression functional
         * @param lambda the lambda, the regression parameter
         */
        DMSystemMatrix(SGPP::base::Grid& grid, SGPP::base::DataMatrix& trainData, SGPP::base::OperationMatrix& C, double lambda);

        /**
         * Std-Destructor
         */
        virtual ~DMSystemMatrix();

        virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

        /**
         * Generates the right hand side of the classification equation
         *
         * @param classes the class information of the training data
         * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
         */
        virtual void generateb(SGPP::base::DataVector& classes, SGPP::base::DataVector& b);
    };

  }
}

#endif /* DMSYSTEMMATRIX_HPP */
