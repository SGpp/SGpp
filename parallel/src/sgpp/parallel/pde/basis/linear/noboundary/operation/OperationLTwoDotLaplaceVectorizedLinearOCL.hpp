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

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    /**
     * Implementation for linear functions of LTwoDotLaplace Operation, linear grids without boundaries using OpenCL
     *
     * @version $HEAD$
     */
    class OperationLTwoDotLaplaceVectorizedLinearOCL: public OperationParabolicPDEMatrixCombined {
      private:
        SGPP::base::GridStorage* storage;
        SGPP::base::DataMatrix* level_;
        SGPP::base::DataMatrix* level_int_;
        SGPP::base::DataMatrix* index_;
        double* lcl_q;
        double* lcl_q_inv;
        SGPP::base::DataVector* lambda;

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
        OperationLTwoDotLaplaceVectorizedLinearOCL(SGPP::base::GridStorage* storage, SGPP::base::DataVector& lambda);

        /**
         * Construtor of OperationLTwoDotLaplaceVectorizedLinearOCL
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationLTwoDotLaplaceVectorizedLinearOCL(SGPP::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~OperationLTwoDotLaplaceVectorizedLinearOCL();

        virtual void mult(SGPP::base::DataVector& alpha,
                          SGPP::base::DataVector& result);

    };

  }

}

#endif /* OPERATIONLTWODOTLAPLACEVECTORIZEDLINEAROCL_HPP */
