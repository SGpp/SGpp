/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Jacob Jepsen (jepsen@diku.dk)

#ifndef OPERATIONLTWODOTLAPLACEVECTORIZEDLINEAROCL_HPP
#define OPERATIONLTWODOTLAPLACEVECTORIZEDLINEAROCL_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/parallel/pde/basis/common/OCLPDEKernels.hpp>
#include <sgpp/parallel/pde/operation/OperationParabolicPDEMatrixCombined.hpp>

namespace sg {
  namespace parallel {

    /**
     * Implementation for linear functions of LTwoDotLaplace Operation, linear grids without boundaries using OpenCL
     *
     * @version $HEAD$
     */
    class OperationLTwoDotLaplaceVectorizedLinearOCL: public OperationParabolicPDEMatrixCombined {
      private:
        sg::base::GridStorage* storage;
        sg::base::DataMatrix* level_;
        sg::base::DataMatrix* level_int_;
        sg::base::DataMatrix* index_;
        double* lcl_q;
        double* lcl_q_inv;
        sg::base::DataVector* lambda;

        OCLPDEKernels OCLPDEKernelsHandle ;
        size_t padding_size;
        size_t sizepad;
        double* subresult;

      public:
        /**
         * Construtor of OperationLTwoDotLaplaceVectorizedLinearOCL
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param lambda Vector which contains pre-factors for every dimension of the operator
         */
        OperationLTwoDotLaplaceVectorizedLinearOCL(sg::base::GridStorage* storage, sg::base::DataVector& lambda);

        /**
         * Construtor of OperationLTwoDotLaplaceVectorizedLinearOCL
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationLTwoDotLaplaceVectorizedLinearOCL(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~OperationLTwoDotLaplaceVectorizedLinearOCL();

        virtual void mult(sg::base::DataVector& alpha,
                          sg::base::DataVector& result);

    };

  }

}

#endif /* OPERATIONLTWODOTLAPLACEVECTORIZEDLINEAROCL_HPP */
