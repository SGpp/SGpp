/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALITERATIVESPCUDALINEAR_HPP
#define OPERATIONMULTIPLEEVALITERATIVESPCUDALINEAR_HPP

#include "base/grid/GridStorage.hpp"
#include "base/tools/SGppStopwatch.hpp"

#include "parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp"

namespace sg {

  namespace parallel {

    /**
     * This class implements OperationMultipleEvalVectorizedSP for a grids with linear basis ansatzfunctions without boundaries
     *
     * However in this case high efficient vector code (CUDA) are generated
     * to implement a iterative OperationB version. In addition cache blocking is used
     * in order to assure a most efficient cache usage.
     *
     * IMPORTANT REMARK:
     * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
     * @li data MUST a have even number of points AND it must be transposed
     * @li result MUST have the same size as data points that should be evaluated
     */
    class OperationMultipleEvalIterativeSPCUDALinear : public sg::parallel::OperationMultipleEvalVectorizedSP {
      public:
        /**
         * Construtor of OperationMultipleEvalVectorizedSP
         *
         * Within the construct DataMatrixSP Level and DataMatrixSP Index are set up.
         * If the grid changes during your calculations and you don't want to create
         * a new instance of this class, you have to call rebuildLevelAndIndex before
         * doing any further mult or multTranspose calls.
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param dataset dataset that should be evaluated
         */
        OperationMultipleEvalIterativeSPCUDALinear(sg::base::GridStorage* storage, sg::base::DataMatrixSP* dataset);

        /**
         * Destructor
         */
        virtual ~OperationMultipleEvalIterativeSPCUDALinear();

        virtual double multVectorized(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result);

        virtual double multTransposeVectorized(sg::base::DataVectorSP& source, sg::base::DataVectorSP& result);

        virtual void rebuildLevelAndIndex(size_t gridFrom = 0, size_t gridTo = std::numeric_limits<size_t>::max());

      protected:
        /// Timer object to handle time measurements
        sg::base::SGppStopwatch* myTimer;
        /// old storage size for deleting buffer on KNF
        size_t old_storage_size;
    };

  }

}

#endif /* OPERATIONMULTIPLEEVALITERATIVESPCUDALINEAR_HPP */
