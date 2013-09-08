/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Jacob Jepsen (jepsen@diku.dk)

#ifndef OPERATIONLTWODOTLAPLACEVECTORIZEDLINEARBOUNDARYOCL_HPP
#define OPERATIONLTWODOTLAPLACEVECTORIZEDLINEARBOUNDARYOCL_HPP

#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/Grid.hpp"
#include "parallel/pde/basis/common/OCLPDEKernels.hpp"
#include "parallel/pde/operation/OperationParabolicPDEMatrixCombined.hpp"

namespace sg {
  namespace parallel {

    /**
     * Implementation for linear functions of LTwoDotLaplace Operation, linear grids with boundaries
     *
     * @version $HEAD$
     */
    class OperationLTwoDotLaplaceVectorizedLinearBoundaryOCL: public OperationParabolicPDEMatrixCombined {
    private:
      sg::base::GridStorage* storage;
      sg::base::DataMatrix* level_;
      sg::base::DataMatrix* level_int_;
      sg::base::DataMatrix* index_;
      double* lcl_q;
      double* lcl_q_inv;
      sg::base::DataVector* lambda;
      OCLPDEKernels OCLPDEKernelsHandle ;

    public:
      /**
       * Construtor of OperationLTwoDotLaplaceLinear
       *
       * @param storage Pointer to the grid's gridstorage obejct
       */
      OperationLTwoDotLaplaceVectorizedLinearBoundaryOCL(sg::base::GridStorage* storage, sg::base::DataVector& lambda);

        /**
         * Construtor of OperationLTwoDotLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param lambda Vector which contains pre-factors for every dimension of the operator
         */
      OperationLTwoDotLaplaceVectorizedLinearBoundaryOCL(sg::base::GridStorage* storage);

      /**
       * Destructor
       */
      virtual ~OperationLTwoDotLaplaceVectorizedLinearBoundaryOCL();

      virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);
    };
  }
}

#endif /* OPERATIONLTWODOTLAPLACEVECTORIZEDOCLLINEARBOUNDARY_HPP */
