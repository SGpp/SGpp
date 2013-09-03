/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Jacob Jepsen (jepsen@diku.dk)

#ifndef OPERATIONLAPLACEVECTORIZEDOCLLINEARBOUNDARY_HPP
#define OPERATIONLAPLACEVECTORIZEDOCLLINEARBOUNDARY_HPP

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
    class OperationLaplaceVectorizedOCLLinearBoundary: public sg::base::OperationMatrix {
    private:
      sg::base::GridStorage* storage;
      sg::base::DataMatrix* level_;
      sg::base::DataMatrix* level_int_;
      sg::base::DataMatrix* index_;
      double* lcl_q;
      double* lcl_q_inv;
      sg::base::DataVector* lambda;
      OCLPDEKernels OCLPDEKernelsHandle;

    public:
      /**
       * Construtor of OperationLaplaceLinear
       *
       * @param storage Pointer to the grid's gridstorage obejct
       */
      OperationLaplaceVectorizedOCLLinearBoundary(sg::base::GridStorage* storage, sg::base::DataVector& lambda);

        /**
         * Construtor of OperationLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param lambda Vector which contains pre-factors for every dimension of the operator
         */
      OperationLaplaceVectorizedOCLLinearBoundary(sg::base::GridStorage* storage);

      /**
       * Destructor
       */
      virtual ~OperationLaplaceVectorizedOCLLinearBoundary();

      virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);
    };

  }

}

#endif /* OPERATIONLAPLACEVECTORIZEDOCLLINEARBOUNDARY_HPP */
