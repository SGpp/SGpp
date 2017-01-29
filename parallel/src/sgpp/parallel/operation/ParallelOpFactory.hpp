// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PARALLEL_OP_FACTORY_HPP
#define PARALLEL_OP_FACTORY_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp>
#include <sgpp/parallel/pde/operation/OperationParabolicPDEMatrixCombined.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>

/*
 * This file contains factory methods for operations.
 */

#include <sgpp/globaldef.hpp>

#include <limits>

namespace sgpp {

namespace op_factory {

/**
 * Factory method, returning an OperationMultipleEvalVectorized for the grid at hand.
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
 * responsible. However, for performance optimizations, implementations might actually calculate
 * (and write to)
 * data that lies beyond these boundaries.
 *
 * @return Pointer to the new OperationMultipleEvalVectorized for the Grid grid
 */
parallel::OperationMultipleEvalVectorized* createOperationMultipleEvalVectorized(
    base::Grid& grid, const parallel::VectorizationType vecType, base::DataMatrix* dataset,
    size_t gridFrom = 0, size_t gridTo = std::numeric_limits<size_t>::max(), size_t datasetFrom = 0,
    size_t datasetTo = std::numeric_limits<size_t>::max());

// #ifdef USE_MPI
/**
 * Factory method, returning an OperationLTwoDotProduct (OperationMatrix) for the grid at hand.
 *
 * This LTwoDotProduct is implemented directly be creating each matrix element on the fly
 *
 * @param grid Grid which is to be used
 * @param vecType selected vectorization
 * @return Pointer to the new OperationMatrix object for the Grid grid
 */
base::OperationMatrix* createOperationLTwoDotProductVectorized(
    base::Grid& grid, const parallel::VectorizationType& vecType);

/**
 * Factory method, returning an OperationLTwoDotLaplace (OperationParabolicPDEMatrixCombined) for
 * the grid at hand.
 * Note: object has to be freed after use.
 *
 * This LTwoDotLaplace operation is implemented directly be creating each matrix element on the fly
 *
 * @param grid Grid which is to be used
 * @param lambda Vector which contains pre-factors for every dimension of the operator
 * @param vecType selected vectorization
 * @return Pointer to the new OperationMatrix object for the Grid grid
 */
parallel::OperationParabolicPDEMatrixCombined*
createOperationLTwoDotLaplaceVectorized(base::Grid& grid, sgpp::base::DataVector& lambda,
                                        const parallel::VectorizationType& vecType);

/**
 * Factory method, returning an OperationLTwoDotLaplace (OperationParabolicPDEMatrixCombined) for
 * the grid at hand.
 *
 * This LTwoDotLaplace operation is implemented directly be creating each matrix element on the fly
 *
 * @param grid Grid which is to be used
 * @param vecType selected vectorization
 * @return Pointer to the new OperationMatrix object for the Grid grid
 */
parallel::OperationParabolicPDEMatrixCombined*
createOperationLTwoDotLaplaceVectorized(base::Grid& grid,
                                        const parallel::VectorizationType& vecType);

/**
 * Factory method, returning an OperationLaplace (OperationMatrix) for the grid at hand.
 *
 * This Laplacian is implemented directly be creating each matrix element on the fly
 *
 * @param grid Grid which is to be used
 * @param lambda Vector which contains pre-factors for every dimension of the operator
 * @param vecType selected vectorization
 * @return Pointer to the new OperationMatrix object for the Grid grid
 */
base::OperationMatrix* createOperationLaplaceVectorized(
    base::Grid& grid, base::DataVector& lambda, const parallel::VectorizationType& vecType);

/**
 * Factory method, returning an OperationLaplace (OperationMatrix) for the grid at hand.
 *
 * This Laplacian is implemented directly be creating each matrix element on the fly
 *
 * @param grid Grid which is to be used
 * @param vecType selected vectorization
 * @return Pointer to the new OperationMatrix object for the Grid grid
 */
base::OperationMatrix* createOperationLaplaceVectorized(
    base::Grid& grid, const parallel::VectorizationType& vecType);
// #endif
}  // namespace op_factory
}  // namespace sgpp

#endif /*PARALLEL_OP_FACTORY_HPP*/
