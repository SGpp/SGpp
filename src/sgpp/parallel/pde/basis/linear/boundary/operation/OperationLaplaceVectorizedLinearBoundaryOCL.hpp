/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Jacob Jepsen (jepsen@diku.dk)

#ifndef OPERATIONLAPLACEVECTORIZEDLINEARBOUNDARYOCL_HPP
#define OPERATIONLAPLACEVECTORIZEDLINEARBOUNDARYOCL_HPP

#include "base/operation/OperationMatrix.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/Grid.hpp"

#include "parallel/pde/basis/common/OCLPDEKernels.hpp"

namespace sg {
  namespace parallel {

    /**
     * Implementation for linear functions of Laplace Operation, linear grids with boundaries
     *
     * @version $HEAD$
     */
    class OperationLaplaceVectorizedLinearBoundaryOCL: public sg::base::OperationMatrix {
      private:
        sg::base::GridStorage* storage;
        sg::base::DataMatrix* level_;
        sg::base::DataMatrix* level_int_;
        sg::base::DataMatrix* index_;
        double* lcl_q;
        double* lcl_q_inv;
        sg::base::DataVector* lambda;
        OCLPDEKernels OCLPDEKernelsHandle;

        void mult_dirichlet(sg::base::DataVector& alpha, sg::base::DataVector& result);

      public:
        /**
         * Construtor of OperationLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param lambda the lambda parameter which is needed in some cases (Black-Scholes) to modify the dimensional local values
         */
        OperationLaplaceVectorizedLinearBoundaryOCL(sg::base::GridStorage* storage, sg::base::DataVector& lambda);

        /**
         * Construtor of OperationLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationLaplaceVectorizedLinearBoundaryOCL(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~OperationLaplaceVectorizedLinearBoundaryOCL();

        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);
    };

  }

}

#endif /* OPERATIONLAPLACEVECTORIZEDLINEARBOUNDARYOCL_HPP */
