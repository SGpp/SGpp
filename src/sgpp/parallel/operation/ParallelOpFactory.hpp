/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)
// @author Alexander Heinecke (alexander.heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef PARALLEL_OP_FACTORY_HPP
#define PARALLEL_OP_FACTORY_HPP

#include "base/grid/Grid.hpp"
#include "base/operation/OperationMatrix.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/datatypes/DataMatrixSP.hpp"

#include "parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp"
#include "parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp"
#include "parallel/tools/TypesParallel.hpp"

/*
 * This file contains factory methods for operations.
 */

namespace sg {

  namespace op_factory {

    /**
     * Factory method, returning an OperationMultipleEvalVectorized for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used for evaluation
     * @param vecType Type of Vectorization used
     * @param dataset Dataset (DataMatrix, data points per row) that is to be evaluated
     * @param datasetFrom lower bound index that describes where to start processing the dataset
     * @param datasetTo upper bound index (exclusive) that describes where to end processing the dataset
     * @param gridFrom lower bound index that describes where to start processing the storage
     * @param gridTo upper bound index (exclusive) that describes where to end processing the storage
     *
     * @note the bounds describe for which part of the dataset/grid the create instance will be
     * responsible. However, for performance optimizations, implementations might actually calculate (and write to)
     * data that lies beyond these boundaries.
     *
     * @return Pointer to the new OperationMultipleEvalVectorized for the Grid grid
     */
    parallel::OperationMultipleEvalVectorized* createOperationMultipleEvalVectorized(base::Grid& grid, const parallel::VectorizationType& vecType, base::DataMatrix* dataset,
        size_t gridFrom = 0, size_t gridTo = std::numeric_limits<size_t>::max(),
        size_t datasetFrom = 0, size_t datasetTo = std::numeric_limits<size_t>::max());
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

    /**
     * Factory method, returning an OperationLaplace (OperationMatrix) for the grid at hand.
     * Note: object has to be freed after use.
     *
     * This Laplacian is implemented directly be creating each matrix element on the fly
     *
     * @param grid Grid which is to be used
     * @return Pointer to the new OperationMatrix object for the Grid grid
     */
    base::OperationMatrix* createOperationLaplaceVectorized(base::Grid& grid);

    /**
     * Factory method, returning an OperationLaplace (OperationMatrix) for the grid at hand.
     * Note: object has to be freed after use.
     *
     * This Laplacian is implemented directly be creating each matrix element on the fly
     *
     * @param grid Grid which is to be used
     * @param lambda Vector which contains pre-factors for every dimension of the operator
     * @return Pointer to the new OperationMatrix object for the Grid grid
     */
    base::OperationMatrix* createOperationLaplaceVectorized(base::Grid& grid, sg::base::DataVector& lambda);
  }

}

#endif /*PARALLEL_OP_FACTORY_HPP*/
