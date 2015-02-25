// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SP_PARALLEL_OP_FACTORY_HPP
#define SP_PARALLEL_OP_FACTORY_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataMatrixSP.hpp>

#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>

/*
 * This file contains factory methods for operations.
 */

#include <sgpp/globaldef.hpp>

// TODO David
#if USE_DOUBLE_PRECISION==0


namespace SGPP {

  namespace op_factory {
    /**
     * Factory method, returning an OperationMultipleEvalVectorizedSP for the grid at hand,
     * single precision.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used for evaluation
     * @param vecType Type of Vectorization used
     * @param dataset Dataset (DataMatrixSP, data points per row) that is to be evaluated
     * @param datasetFrom lower bound index that describes where to start processing the dataset
     * @param datasetTo upper bound index (exclusive) that describes where to end processing the dataset
     * @param gridFrom lower bound index that describes where to start processing the storage
     * @param gridTo upper bound index (exclusive) that describes where to end processing the storage
     *
     * @note the bounds describe for which part of the dataset/grid the create instance will be
     * responsible. However, for performance optimizations, implementations might actually calculate (and write to)
     * data that lies beyond these boundaries.
     *
     * @return Pointer to the new OperationMultipleEvalVectorizedSP for the Grid grid
     */
    parallel::OperationMultipleEvalVectorizedSP* createOperationMultipleEvalVectorizedSP(base::Grid& grid, const parallel::VectorizationType& vecType, base::DataMatrixSP* dataset,
        size_t gridFrom = 0, size_t gridTo = std::numeric_limits<size_t>::max(),
        size_t datasetFrom = 0, size_t datasetTo = std::numeric_limits<size_t>::max());
  }

}

#endif

#endif /*SP_PARALLEL_OP_FACTORY_HPP*/
