// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALVECTORIZEDSP_HPP
#define OPERATIONMULTIPLEEVALVECTORIZEDSP_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/parallel/datadriven/tools/LevelIndexMaskOffsetHelper.hpp>

#include <sgpp/globaldef.hpp>

#include <limits>

namespace sgpp {
namespace parallel {

/**
 * @brief Interface for multiplication with Matrices @f$B@f$ and @f$B^T@f$.
 *
 * If there are @f$N@f$ basis functions, @f$\{\varphi(\vec{x})\}_{i=1,\ldots,N}@f$ and @f$m@f$ data
 * points
 *
 * This class defines an interface similar to OperationMultipleEval in order to support SIMD
 * architectures
 * for datadriven task (multiple function evaluations, classification, regression). Target
 * architectures may be Intel SSE, Intel AVX, nVidia CUDA, OpenCL.
 */
class OperationMultipleEvalVectorizedSP {
 protected:
  /// Pointer to the grid's GridStorage object
  sgpp::base::GridStorage& storage_;
  /// Pointer to the dataset that should be evaluated on the grid
  sgpp::base::DataMatrixSP* dataset_;
  /// Member to store the sparse grid's levels for better vectorization
  sgpp::base::DataMatrixSP* level_;
  /// Member to store the sparse grid's indices for better vectorization
  sgpp::base::DataMatrixSP* index_;
  /// Member to store for masks per grid point for better vectorization of modlinear operations
  sgpp::base::DataMatrixSP* mask_;
  /// Member to store offsets per grid point for better vecotrization of modlinear operations
  sgpp::base::DataMatrixSP* offset_;
  /// Timer object to handle time measurements
  sgpp::base::SGppStopwatch* myTimer_;

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
  OperationMultipleEvalVectorizedSP(sgpp::base::GridStorage* storage,
                                    sgpp::base::DataMatrixSP* dataset);

  /**
   * Destructor
   *
   * cleans up level_ and index_ members
   */
  virtual ~OperationMultipleEvalVectorizedSP();

  /**
   * Multiplication of @f$B^T@f$ with vector @f$\alpha@f$
   *
   * IMPORTANT REMARK:
   * In order to use this routine you have to keep following points in mind (for multVectorized and
   * multTransposeVectorized):
   * @li data MUST a have even number of points AND it must be transposed
   * @li result MUST have the same size as data points that should be evaluated
   *
   * @param alpha vector, to which @f$B@f$ is applied. Typically the coefficient vector
   * @param result the result vector of the matrix vector multiplication
   */
  virtual double multVectorized(sgpp::base::DataVectorSP& alpha,
                                sgpp::base::DataVectorSP& result) = 0;

  /**
   * Multiplication of @f$B@f$ with vector @f$\alpha@f$
   *
   * IMPORTANT REMARK:
   * In order to use this routine you have to keep following points in mind (for multVectorized and
   * multTransposeVectorized):
   * @li data MUST a have even number of points AND it must be transposed
   * @li result MUST have the same size as data points that should be evaluated
   *
   * @param source vector, to which @f$B^T@f$ is applied. Typically the coefficient vector
   * @param result the result vector of the matrix vector multiplication
   */
  virtual double multTransposeVectorized(sgpp::base::DataVectorSP& source,
                                         sgpp::base::DataVectorSP& result) = 0;

  /**
   * rebuilds the DataMatrix for Level and Index in Derivatives
   * needed for vectorization.
   *
   * @param gridFrom  the new lower bound for the grid this process should handle (inclusive)
   * @param gridTo  the new upper bound for the grid this process should handle (exclusive)
   *
   * if you don't supply the optional parameters, this process will handle the complete grid
   */
  virtual void rebuildLevelAndIndex(size_t gridFrom = 0,
                                    size_t gridTo = std::numeric_limits<size_t>::max()) = 0;

  // we have to make sure that all implemented Kernel Types can see the fields they require.
  // this has to be done here and only here because friends are checked at compile time and we use
  // a pointer of this class to access the respective objects
  friend struct LevelIndexMaskOffsetHelperSP::rebuild<Standard, OperationMultipleEvalVectorizedSP>;
  friend struct LevelIndexMaskOffsetHelperSP::rebuild<Mask, OperationMultipleEvalVectorizedSP>;
};
}  // namespace parallel
}  // namespace sgpp

#endif /* OPERATIONMULTIPLEEVALVECTORIZEDSP_HPP */
