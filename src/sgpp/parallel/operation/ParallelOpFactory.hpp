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
#include "base/datatypes/DataMatrix.hpp"
#include "base/datatypes/DataMatrixSP.hpp"

#include "parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp"
#include "parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp"
#include "parallel/tools/TypesParallel.hpp"

/*
 * This file contains factory methods for operations.
 */
 
namespace sg
{

namespace op_factory
{

  /**
   * Factory method, returning an OperationMultipleEvalVectorized for the grid at hand.
   * Note: object has to be freed after use.
   *
   * @param grid Grid which is to be used for evaluation
   * @param vecType Type of Vectorization used. Currently supported: SSE, AVX
   * @param dataset Dataset (DataMatrix, data points per row) that is to be evaluated
   * @param multFrom lower bound index that describes where to start processing the dataset
   * @param multTo upper bound index (exclusive) that describes where to end processing the dataset
   * @param multTransposeFrom lower bound index that describes where to start processing the storage
   * @param multTransposeTo upper bound index (exclusive) that describes where to end processing the storage
   *
   * @note the bounds describe for which part of the dataset/grid the create instance will be
   * responsible. However, for performance optimizations, implementations might actually calculate (and write to)
   * data that lies beyond these boundaries.
   *
   * @return Pointer to the new OperationMultipleEvalVectorized for the Grid grid
   */
  parallel::OperationMultipleEvalVectorized* createOperationMultipleEvalVectorized(base::Grid& grid, const parallel::VectorizationType& vecType, base::DataMatrix* dataset,
                                                                                   int gridFrom = 0, int gridTo = -1, int datasetFrom = 0, int datasetTo = -1);
  /**
   * Factory method, returning an OperationMultipleEvalVectorizedSP for the grid at hand,
   * single precision.
   * Note: object has to be freed after use.
   *
   * @param grid Grid which is to be used for evaluation
   * @param vecType Type of Vectorization used. Currently supported: SSE, AVX
   * @param dataset Dataset (DataMatrixSP, data points per row) that is to be evaluated
   * @param multFrom lower bound index that describes where to start processing the dataset
   * @param multTo upper bound index (exclusive) that describes where to end processing the dataset
   * @param multTransposeFrom lower bound index that describes where to start processing the storage
   * @param multTransposeTo upper bound index (exclusive) that describes where to end processing the storage
   *
   * @note the bounds describe for which part of the dataset/grid the create instance will be
   * responsible. However, for performance optimizations, implementations might actually calculate (and write to)
   * data that lies beyond these boundaries.
   *
   * @return Pointer to the new OperationMultipleEvalVectorizedSP for the Grid grid
   */
  parallel::OperationMultipleEvalVectorizedSP* createOperationMultipleEvalVectorizedSP(base::Grid& grid, const parallel::VectorizationType& vecType, base::DataMatrixSP* dataset,
																					   int gridFrom = 0, int gridTo = -1, int datasetFrom = 0, int datasetTo = -1);
}

}

#endif /*PARALLEL_OP_FACTORY_HPP*/
