/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Jacob Jepsen (jepsen@diku.dk)

#ifndef OPERATIONLTWODOTLAPLACEVECTORIZEDOCLLINEARBOUNDARY_HPP
#define OPERATIONLTWODOTLAPLACEVECTORIZEDOCLLINEARBOUNDARY_HPP

#include "base/operation/OperationMatrix.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/Grid.hpp"
#include "parallel/pde/basis/common/OCLPDEKernels.hpp" 

namespace sg {
  namespace parallel {

    /**
     * Implementation for linear functions of LTwoDotLaplace Operation, linear grids with boundaries
     *
     * @version $HEAD$
     */
    class OperationLTwoDotLaplaceVectorizedOCLLinearBoundary: public sg::base::OperationMatrix {
    private:
      sg::base::GridStorage* storage;
      sg::base::DataMatrix* level_;
      sg::base::DataMatrix* level_int_;
      sg::base::DataMatrix* index_;
      double* lcl_q;
      double* lcl_q_inv;
      sg::base::DataVector* lambda;
      double TimestepCoeff;

      OCLPDEKernels OCLPDEKernelsHandle ;

    public:
      /**
       * Construtor of OperationLTwoDotLaplaceLinear
       *
       * @param storage Pointer to the grid's gridstorage obejct
       */
      OperationLTwoDotLaplaceVectorizedOCLLinearBoundary(sg::base::GridStorage* storage, sg::base::DataVector& lambda);

        /**
         * Construtor of OperationLTwoDotLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param lambda Vector which contains pre-factors for every dimension of the operator
         */
      OperationLTwoDotLaplaceVectorizedOCLLinearBoundary(sg::base::GridStorage* storage);

      /**
       * Destructor
       */
      virtual ~OperationLTwoDotLaplaceVectorizedOCLLinearBoundary();

      virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);
      
      /**
       * set timestep coefficient
       *
       * @param newTimestepCoeff new timestep coefficient
       */
      void setTimestepCoeff(double newTimestepCoeff);
      
      /**
       * get current timestep coefficient
       *
       * @return current timestep coefficient
       */
      double getTimestepCoeff();
    };
  }
}

#endif /* OPERATIONLTWODOTLAPLACEVECTORIZEDOCLLINEARBOUNDARY_HPP */
