/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALITERATIVEARBBLINEAR_HPP
#define OPERATIONMULTIPLEEVALITERATIVEARBBLINEAR_HPP

#include "parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp"
#include "parallel/datadriven/basis/common/ArBBKernels.hpp"
#include "parallel/datadriven/basis/common/ArBBKernels2D.hpp"
#include "parallel/datadriven/basis/common/ArBBKernels4D.hpp"
#include "parallel/datadriven/basis/common/ArBBKernels5D.hpp"
#include "parallel/datadriven/basis/common/ArBBKernels10D.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/tools/SGppStopwatch.hpp"

namespace sg {
  namespace parallel {

    /**
     * This class implements sg::base::OperationMultipleEval for a grids with linear basis ansatzfunctions without boundaries
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
    class OperationMultipleEvalIterativeArBBLinear : public sg::parallel::OperationMultipleEvalVectorized {
      public:
        /**
         * Construtor of sg::base::OperationMultipleEvalLinear
         *
         * Within the construct sg::base::DataMatrix Level and sg::base::DataMatrix Index are set up.
         * If the grid changes during your calculations and you don't want to create
         * a new instance of this class, you have to call rebuildLevelAndIndex before
         * doing any further mult or multTranspose calls.
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param dataset dataset that should be evaluated on the grid
         */
        OperationMultipleEvalIterativeArBBLinear(sg::base::GridStorage* storage, sg::base::DataMatrix* dataset);

        /**
         * Destructor
         */
        virtual ~OperationMultipleEvalIterativeArBBLinear();

        virtual double multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result);

        virtual double multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result);

        virtual void rebuildLevelAndIndex();

      protected:
        /// Pointer to the grid's gridstorage object
        sg::base::GridStorage* storage;
        /// Timer object to handle time measurements
        sg::base::SGppStopwatch* myTimer;
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

#endif /* OPERATIONMULTIPLEEVALITERATIVEARBBLINEAR_HPP */
