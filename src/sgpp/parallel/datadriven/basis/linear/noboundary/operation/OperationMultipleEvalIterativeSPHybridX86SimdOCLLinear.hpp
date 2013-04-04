/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALITERATIVESPHYBRIDX86SIMDOCLLINEAR_HPP
#define OPERATIONMULTIPLEEVALITERATIVESPHYBRIDX86SIMDOCLLINEAR_HPP

#include "parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp"
#include "parallel/datadriven/basis/common/OCLKernels.hpp"
#include "parallel/tools/TwoPartitionAutoTuning.hpp"

namespace sg {
  namespace parallel {

    /**
     * This class implements OperationMultipleEvalSP for a grids with linear basis ansatzfunctions without boundaries.
     *
     * However in this case highly efficient vector code (OpenCL and SSE/AVX intrinsics) is generated
     * to implement a iterative OperationB version. In addition cache blocking is used
     * in order to assure a most efficient cache usage.
     *
     * IMPORTANT REMARK:
     * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
     * @li data MUST a have even number of points AND it must be transposed
     * @li result MUST have the same size as data points that should be evaluated
     */
    class OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear : public sg::parallel::OperationMultipleEvalVectorizedSP {
      public:
        /**
         * Construtor of OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear
         *
         * Within the construct sg::base::DataMatrixSP Level and sg::base::DataMatrixSP Index are set up.
         * If the grid changes during your calculations and you don't want to create
         * a new instance of this class, you have to call rebuildLevelAndIndex before
         * doing any further mult or multTranspose calls.
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param dataset dataset that should be evaluated
         */
        OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear(sg::base::GridStorage* storage, sg::base::DataMatrixSP* dataset);

        /**
         * Destructor
         */
        virtual ~OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear();

        virtual double multVectorized(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result);

        virtual double multTransposeVectorized(sg::base::DataVectorSP& source, sg::base::DataVectorSP& result);

        virtual void rebuildLevelAndIndex();

      protected:
        /// Object to access the OCL Kernel
        OCLKernels* myOCLKernels;
        /// Autotuning object for mult routine
        sg::parallel::TwoPartitionAutoTuning* _tuningMult;
        /// Autotuning object for mult trans routine
        sg::parallel::TwoPartitionAutoTuning* _tuningMultTrans;
    };

  }
}

#endif /* OPERATIONMULTIPLEEVALITERATIVESPHYBRIDX86SIMDOCLLINEAR_HPP */
