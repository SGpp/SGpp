// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALITERATIVESPARBBLINEAR_HPP
#define OPERATIONMULTIPLEEVALITERATIVESPARBBLINEAR_HPP

#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp>
#include <sgpp/parallel/datadriven/basis/common/ArBBKernels.hpp>
#include <sgpp/parallel/datadriven/basis/common/ArBBKernels2D.hpp>
#include <sgpp/parallel/datadriven/basis/common/ArBBKernels4D.hpp>
#include <sgpp/parallel/datadriven/basis/common/ArBBKernels5D.hpp>
#include <sgpp/parallel/datadriven/basis/common/ArBBKernels10D.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    /**
     * This class implements OperationMultipleEvalSP for a grids with linear basis ansatzfunctions without boundaries
     *
     * However in this case high efficient vector code (Intel Array Building Blocks) is generated
     * to implement a iterative OperationB version. In addition cache blocking is used
     * in order to assure a most efficient cache usage.
     *
     * IMPORTANT REMARK:
     * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
     * @li data MUST a have even number of points AND it must be transposed
     * @li result MUST have the same size as data points that should be evaluated
     */
    class OperationMultipleEvalIterativeSPArBBLinear : public SGPP::parallel::OperationMultipleEvalVectorizedSP {
      public:
        /**
         * Construtor of OperationBLinearSP
         *
         * Within the construct SGPP::base::DataMatrixSP Level and SGPP::base::DataMatrixSP Index are set up.
         * If the grid changes during your calculations and you don't want to create
         * a new instance of this class, you have to call rebuildLevelAndIndex before
         * doing any further mult or multTranspose calls.
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param dataset dataset that should be evaluated
         */
        OperationMultipleEvalIterativeSPArBBLinear(SGPP::base::GridStorage* storage, SGPP::base::DataMatrixSP* dataset);

        /**
         * Destructor
         */
        virtual ~OperationMultipleEvalIterativeSPArBBLinear();

        virtual double multVectorized(SGPP::base::DataVectorSP& alpha, SGPP::base::DataVectorSP& result);

        virtual double multTransposeVectorized(SGPP::base::DataVectorSP& source, SGPP::base::DataVectorSP& result);

        virtual void rebuildLevelAndIndex();

      protected:
        /// Pointer to the grid's gridstorage object
        SGPP::base::GridStorage* storage;
        /// Timer object to handle time measurements
        SGPP::base::SGppStopwatch* myTimer;
        /// Object to access the OCL Kernel
        ArBBKernels* myArBBKernels;
        /// Object to access the ArBB 2D Kernel
        ArBBKernels2D* myArBBKernels2D;
        /// Object to access the ArBB 4D Kernel
        ArBBKernels4D* myArBBKernels4D;
        /// Object to access the ArBB 5D Kernel
        ArBBKernels5D* myArBBKernels5D;
        /// Object to access the ArBB 10D Kernel
        ArBBKernels10D* myArBBKernels10D;
    };

  }
}

#endif /* OPERATIONMULTIPLEEVALITERATIVESPARBBLINEAR_HPP */