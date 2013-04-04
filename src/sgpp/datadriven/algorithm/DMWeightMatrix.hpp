/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Zhongwen Song (songz@in.tum.de)
// @author Benjamin (pehersto@in.tum.de)

#ifndef DMWEIGHTMATRIX_HPP
#define DMWEIGHTMATRIX_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/grid/Grid.hpp"
#include "base/operation/OperationMultipleEval.hpp"
#include "base/operation/OperationMatrix.hpp"

namespace sg {
  namespace datadriven {

    /**
     * Class that implements the virtual class OperationMatrix for the
     * application of classification for the Systemmatrix with weight
     */
    class DMWeightMatrix : public sg::base::OperationMatrix {
      private:
        /// the lambda, the regularisation parameter
        double lamb;
        /// sg::base::OperationMatrix, the regularisation mehtod
        sg::base::OperationMatrix* C;
        /// OperationB for calculating the data matrix
        sg::base::OperationMultipleEval* B;
        /// Pointer to the data vector
        sg::base::DataMatrix* data;
        /// Pointer to the weight vector
        sg::base::DataVector* weight;

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param trainData reference to sg::base::DataVector that contains the training data
         * @param C the regression functional
         * @param lambda the lambda, the regression parameter
         * @param w the weights to the training data
         */
        DMWeightMatrix(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, sg::base::OperationMatrix& C, double lambda, sg::base::DataVector& w);

        /**
         * Std-Destructor
         */
        virtual ~DMWeightMatrix();

        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

        /**
         * Generates the right hand side of the classification equation
         *
         * @param classes the class information of the training data
         * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
         */
        void generateb(sg::base::DataVector& classes, sg::base::DataVector& b);
    };

  }
}

#endif /* DMWEIGHTMATRIX_HPP */
