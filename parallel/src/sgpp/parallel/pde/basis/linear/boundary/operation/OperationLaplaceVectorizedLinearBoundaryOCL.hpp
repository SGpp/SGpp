// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONLAPLACEVECTORIZEDLINEARBOUNDARYOCL_HPP
#define OPERATIONLAPLACEVECTORIZEDLINEARBOUNDARYOCL_HPP

#include <sgpp/base/operation/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/parallel/pde/basis/common/OCLPDEKernels.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    /**
     * Implementation for linear functions of Laplace Operation, linear grids with boundaries
     *
     * @version $HEAD$
     */
    class OperationLaplaceVectorizedLinearBoundaryOCL: public SGPP::base::OperationMatrix {
      private:
        SGPP::base::GridStorage* storage;
        SGPP::base::DataMatrix* level_;
        SGPP::base::DataMatrix* level_int_;
        SGPP::base::DataMatrix* index_;
        double* lcl_q;
        double* lcl_q_inv;
        SGPP::base::DataVector* lambda;
        OCLPDEKernels OCLPDEKernelsHandle;

        void mult_dirichlet(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

      public:
        /**
         * Construtor of OperationLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param lambda the lambda parameter which is needed in some cases (Black-Scholes) to modify the dimensional local values
         */
        OperationLaplaceVectorizedLinearBoundaryOCL(SGPP::base::GridStorage* storage, SGPP::base::DataVector& lambda);

        /**
         * Construtor of OperationLaplaceLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationLaplaceVectorizedLinearBoundaryOCL(SGPP::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~OperationLaplaceVectorizedLinearBoundaryOCL();

        virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);
    };

  }

}

#endif /* OPERATIONLAPLACEVECTORIZEDLINEARBOUNDARYOCL_HPP */