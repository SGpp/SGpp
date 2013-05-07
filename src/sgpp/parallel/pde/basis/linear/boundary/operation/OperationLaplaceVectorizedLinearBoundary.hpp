/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONLAPLACEVECTORIZEDLINEARBOUNDARY_HPP
#define OPERATIONLAPLACEVECTORIZEDLINEARBOUNDARY_HPP

#include "base/operation/OperationMatrix.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/Grid.hpp"

namespace sg {
  namespace parallel {

    /**
     * Implementation for linear functions of Laplace Operation, linear grids with boundaries
     *
     * @version $HEAD$
     */
    class OperationLaplaceVectorizedLinearBoundary: public sg::base::OperationMatrix {
      private:
        sg::base::GridStorage* storage;
        sg::base::DataMatrix* level_;
        sg::base::DataMatrix* level_int_;
        sg::base::DataMatrix* index_;
        double* lcl_q;
        double* lcl_q_inv;
        sg::base::DataVector* lambda_;

        double gradient_dirichlet(size_t i, size_t j, size_t dim);

        double l2dot_dirichlet(size_t i, size_t j, size_t dim);

        void mult_dirichlet(sg::base::DataVector& alpha, sg::base::DataVector& result);

      public:
        /**
         * Construtor of OperationLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationLaplaceVectorizedLinearBoundary(sg::base::GridStorage* storage);

        /**
         * Construtor of OperationLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param lambda Vector which contains pre-factors for every dimension of the operator
         */
        OperationLaplaceVectorizedLinearBoundary(sg::base::GridStorage* storage, sg::base::DataVector& lambda);

        /**
         * Destructor
         */
        virtual ~OperationLaplaceVectorizedLinearBoundary();

        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);
    };

  }

}

#endif /* OPERATIONLAPLACEVECTORIZEDLINEARBOUNDARY_HPP */
