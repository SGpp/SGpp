// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef __AVX__

#pragma once

#include <iostream>
#include <map>
#include <vector>

//////////////////////////////////////////////////////////////////////
// Caution: Subspace-skipping is disabled by default for this kernel
//////////////////////////////////////////////////////////////////////

#include <omp.h>

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/AbstractOperationMultipleEvalSubspace.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimpleParameters.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/SubspaceNodeSimple.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Multiple evaluation operation that uses the subspace structure to save work
 * compared to the naive or streaming variants.
 * Verification and debugging variant, should not be used. Instead, use the
 * @ref OperationMultipleEvalSubspaceCombined.
 */
class OperationMultipleEvalSubspaceSimple : public AbstractOperationMultipleEvalSubspace {
 private:
  size_t dim = -1;
  size_t maxLevel = 0;

  size_t* allSubspaces = NULL;
  size_t subspaceCount = -1;
  size_t subspaceSize = -1;

  size_t totalGridPoints = 0;
  double* allSurplusses = nullptr;
  std::map<uint32_t, uint32_t> allSurplussesIndexMap;

  void prepareSubspaceIterator();

  void createFlatStorage();

  void setSurplus(std::vector<size_t>& level, std::vector<size_t>& maxIndices,
                  std::vector<size_t>& index, double value);

  void getSurplus(std::vector<size_t>& level, std::vector<size_t>& maxIndices,
                  std::vector<size_t>& index, double& value, bool& isVirtual);

  void setCoefficients(base::DataVector& surplusVector);

  void unflatten(base::DataVector& result);

  // NOLINT
  size_t flattenIndex(size_t dim, std::vector<size_t>& maxIndices, std::vector<size_t>& index);

  // NOLINT
  size_t flattenIndex(size_t* intermediates, size_t dim, size_t* maxIndicesPtr, size_t* indexPtr,
                      size_t toRecalc);

  size_t flattenLevel(size_t dim, size_t maxLevel, std::vector<size_t>& level);

  static inline size_t calculateIndexComponent(double unadjusted) {
    // implies flooring
    size_t rounded = static_cast<size_t>(unadjusted);

    size_t mask = 0x1;
    size_t sign = mask ^ (mask & rounded);

    size_t componentIndex = rounded + sign;
    return componentIndex;
  }

 public:
  /**
   * Creates a new instance of the OperationMultipleEvalSubspaceSimple class.
   *
   * @param grid grid to be evaluated
   * @param dataset set of evaluation points
   */
  OperationMultipleEvalSubspaceSimple(base::Grid& grid, base::DataMatrix& dataset);

  /**
   * Destructor
   */
  ~OperationMultipleEvalSubspaceSimple();

  /**
   * Updates the internal data structures to reflect changes to the grid, e.g. due to refinement.
   *
   */
  void prepare() override;

  /**
   * Internal eval operator, should not be called directly.
   *
   * @see OperationMultipleEval
   *
   * @param alpha surplusses of the grid
   * @param result will contain the evaluation results for the given range.
   * @param start_index_data beginning of the range to evaluate
   * @param end_index_data end of the range to evaluate
   */
  void multTransposeImpl(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                         const size_t start_index_data, const size_t end_index_data) override;

  /**
   * Internal mult operator, should not be called directly.
   *
   * @see OperationMultipleEval
   *
   * @param source source operand for the operator
   * @param result stores the result
   * @param start_index_data beginning of the range to process
   * @param end_index_data end of the range to process
   */
  void multImpl(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                const size_t start_index_data, const size_t end_index_data) override;

  /**
   * Alignment required by the vector instruction set SG++ is compiled with.
   *
   * @result alignment requirement
   */
  size_t getAlignment() override;

  /**
   * Name of the implementation, useful for benchmarking different implementation approaches.
   *
   * @result name of the implementation
   */
  std::string getImplementationName() override;
};
}
}

#endif
