/* ****************************************************************************
* Copyright (C) 2011-2013 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef OPERATIONMULTIPLEEVALVECTORIZED_HPP
#define OPERATIONMULTIPLEEVALVECTORIZED_HPP

#include "base/grid/GridStorage.hpp"
#include "base/operation/OperationMatrix.hpp"
#include "base/tools/SGppStopwatch.hpp"
#include "parallel/datadriven/tools/LevelIndexMaskOffsetHelper.hpp"

namespace sg {
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
        sg::base::GridStorage* storage_;
        /// Pointer to the dataset that should be evaluated on the grid
        sg::base::DataMatrix* dataset_;
        /// Member to store the sparse grid's levels for better vectorization
        sg::base::DataMatrix* level_;
        /// Member to store the sparse grid's indices for better vectorization
        sg::base::DataMatrix* index_;
        /// Member to store for masks per grid point for better vectorization of modlinear operations
        sg::base::DataMatrix* mask_;
        /// Member to store offsets per grid point for better vecotrization of modlinear operations
        sg::base::DataMatrix* offset_;
        /// Timer object to handle time measurements
        sg::base::SGppStopwatch* myTimer_;

        int m_gridFrom;
        int m_gridTo;
        int m_datasetFrom;
        int m_datasetTo;

      public:
        /**
         * Constructor
         *
         * @param storage GridStorage object used in this operation
         * @param dataset data set that should be evaluated on the sparse grid
         */
        OperationMultipleEvalVectorized(sg::base::GridStorage* storage, sg::base::DataMatrix* dataset);

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
        virtual double multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;

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
        virtual double multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result) = 0;

        /**
         * rebuilds the DataMatrix for Level and Index in Derivatives
         * needed for vectorization.
         */
        virtual void rebuildLevelAndIndex() = 0;

        /**
         * @brief updates the compute boundaries for the grid, after this has been resized
         *
         * @todo for now, the default implementation does nothing. perhaps remove default implementation.
         * it would be an idea to integrate this with rebuildLevelAndIndex
         *
         * @param gridFrom
         * @param gridTo
         */
        virtual void updateGridComputeBoundaries(int gridFrom, int gridTo) {}

        // we have to make sure that all implemented Kernel Types can see the fields they require.
        // this has to be done here and only here because friends are checked at compile time and we use
        // a pointer of this class to access the respective objects
        friend struct LevelIndexMaskOffsetHelper::rebuild<Standard, OperationMultipleEvalVectorized>;
        friend struct LevelIndexMaskOffsetHelper::rebuild<Mask, OperationMultipleEvalVectorized>;
    };

  }
}

#endif /* OPERATIONMULTIPLEEVALVECTORIZED_HPP */
