/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Fabian Franzelin <franzeli@in.tum.de>,
//         Benjamin Peherstorfer <pehersto@in.tum.de>

#ifndef DENSITYSYSTEMMATRIX_HPP
#define DENSITYSYSTEMMATRIX_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/OperationMultipleEval.hpp>
#include <sgpp/base/operation/OperationMatrix.hpp>

namespace sg {
  namespace datadriven {

    /**
     * Class that implements the virtual class OperationMatrix for the
     * application of classification for the Systemmatrix by using a
     * density function
     */
    class DensitySystemMatrix : public sg::base::OperationMatrix {
      private:
        /// the lambda, the regularisation parameter
        double lambda;
        /// Operation A for calculating the data matrix
        /// (L2 Dot-Product of basis functions)
        sg::base::OperationMatrix* A;
        /// OperationB for calculating the data matrix
        sg::base::OperationMultipleEval* B;
        /// OperationMatrix, the regularisation method
        sg::base::OperationMatrix* C;
        /// Training data
        sg::base::DataMatrix* data;

      public:
        /**
         * Std-Constructor
         *
         * @param grid  reference to the sparse grid
         * @param trainData reference to DataVector that contains the training data
         * @param C the regression functional
         * @param lambda the regression parameter
         */
        DensitySystemMatrix(sg::base::Grid& grid, sg::base::DataMatrix& trainData,
                            sg::base::OperationMatrix& C, double lambda);

        /**
         * Generates the left hand side of the classification equation
         *
         * @param alpha parameters for the sparse grid functions
         * @param result reference to the vector which will contain the result
         */
        void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

        /**
         * Generates the right hand side of the classification equation
         *
         * @param b reference to the vector which will contain the result of the
         * matrix vector multiplication on the rhs
         */
        void generateb(sg::base::DataVector& b);

        /**
         * Std-Destructor
         */
        virtual ~DensitySystemMatrix();
    };

  }
}

#endif /* DENSITYSYSTEMMATRIX_HPP */
