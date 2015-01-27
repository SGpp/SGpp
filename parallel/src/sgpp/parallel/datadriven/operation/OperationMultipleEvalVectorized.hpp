// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALVECTORIZED_HPP
#define OPERATIONMULTIPLEEVALVECTORIZED_HPP

#include <limits>

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/parallel/datadriven/tools/LevelIndexMaskOffsetHelper.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    /**
     * @brief Interface for multiplication with Matrices @f$B@f$ and @f$B^T@f$.
     *
     * If there are @f$N@f$ basis functions, @f$\{\varphi(\vec{x})\}_{i=1,\ldots,N}@f$ and @f$m@f$ data points
     *
     * This class defines an interface similar to OperationMultipleEval in order to support SIMD architectures
     * for datadriven task (multiple function evaluations, classification, regression). Target
     * architectures may be Intel SSE, Intel AVX, nVidia CUDA, OpenCL.
     */
    class OperationMultipleEvalVectorized {
      protected:
        /// Pointer to the grid's gridstorage object
        SGPP::base::GridStorage* storage_;
        /// Pointer to the dataset that should be evaluated on the grid
        SGPP::base::DataMatrix* dataset_;
        /// Member to store the sparse grid's levels for better vectorization
        SGPP::base::DataMatrix* level_;
        /// Member to store the sparse grid's indices for better vectorization
        SGPP::base::DataMatrix* index_;
        /// Member to store for masks per grid point for better vectorization of modlinear operations
        SGPP::base::DataMatrix* mask_;
        /// Member to store offsets per grid point for better vecotrization of modlinear operations
        SGPP::base::DataMatrix* offset_;
        /// Timer object to handle time measurements
        SGPP::base::SGppStopwatch* myTimer_;

        /// @todo add boundaries to constructor of this class
        size_t m_gridFrom;
        size_t m_gridTo;
        size_t m_datasetFrom;
        size_t m_datasetTo;

      public:
        /**
         * Constructor
         *
         * @param storage GridStorage object used in this operation
         * @param dataset data set that should be evaluated on the sparse grid
         */
        OperationMultipleEvalVectorized(SGPP::base::GridStorage* storage, SGPP::base::DataMatrix* dataset);

        /**
         * Destructor
         *
         * cleans up level_ and index_ members
         */
        virtual ~OperationMultipleEvalVectorized();

        /**
         * Multiplication of @f$B^T@f$ with vector @f$\alpha@f$
         *
         * IMPORTANT REMARK:
         * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
         * @li data MUST a have even number of points AND it must be transposed
         * @li result MUST have the same size as data points that should be evaluated
         *
         * @param alpha vector, to which @f$B@f$ is applied. Typically the coefficient vector
         * @param result the result vector of the matrix vector multiplication
         */
        virtual double multVectorized(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;

        /**
         * Multiplication of @f$B@f$ with vector @f$\alpha@f$
         *
         * IMPORTANT REMARK:
         * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
         * @li data MUST a have even number of points AND it must be transposed
         * @li result MUST have the same size as data points that should be evaluated
         *
         * @param source vector, to which @f$B^T@f$ is applied. Typically the coefficient vector
         * @param result the result vector of the matrix vector multiplication
         */
        virtual double multTransposeVectorized(SGPP::base::DataVector& source, SGPP::base::DataVector& result) = 0;

        /**
         * rebuilds the DataMatrix for Level and Index in Derivatives
         * needed for vectorization.
         *
         * @param gridFrom  the new lower bound for the grid this process should handle (inclusive)
         * @param gridTo  the new upper bound for the grid this process should handle (exclusive)
         *
         * if you don't supply the optional parameters, this process will handle the complete grid
         */
        virtual void rebuildLevelAndIndex(size_t gridFrom = 0, size_t gridTo = std::numeric_limits<size_t>::max()) = 0;

        // we have to make sure that all implemented Kernel Types can see the fields they require.
        // this has to be done here and only here because friends are checked at compile time and we use
        // a pointer of this class to access the respective objects
        friend struct LevelIndexMaskOffsetHelper::rebuild<Standard, OperationMultipleEvalVectorized>;
        friend struct LevelIndexMaskOffsetHelper::rebuild<Mask, OperationMultipleEvalVectorized>;
    };

  }
}

#endif /* OPERATIONMULTIPLEEVALVECTORIZED_HPP */