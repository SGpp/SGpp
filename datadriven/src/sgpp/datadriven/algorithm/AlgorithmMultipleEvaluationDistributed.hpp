/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * AlgorithmMultipleEvaluationDistributed.hpp
 *
 * Created on: Mar 23, 2019
 *     Author: Jan Schopohl
 */

#pragma once

#include <sgpp/base/algorithm/AlgorithmEvaluation.hpp>
#include <sgpp/base/algorithm/AlgorithmEvaluationTransposed.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>
#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

/**
 * Abstract implementation for distributed multiple function evaluations.
 *
 * Based on AlgorithmMultipleEvaluation.hpp
 *
 * In Data Mining to operators are needed: mass evaluation and transposed evaluation, referenced
 * in literature as matrix vector products with matrices B^T (mass evaluation) and B
 * (transposed evaluation).
 *
 * If there are @f$N@f$ basis functions @f$\varphi(\vec{x})@f$ and @f$m@f$ data points, then B is a
 * (mxN) matrix, with
 * @f[ (B)_{j,i} = \varphi_i(x_j). @f]
 *
 */
template <class BASIS>
class AlgorithmMultipleEvaluationDistributed {
 public:
  /**
   * Performs a distributed transposed mass evaluation
   *
   * @param storage GridStorage object that contains the grid's points information
   * @param basis a reference to a class that implements a specific basis
   * @param source the coefficients of the grid points
   * @param x the d-dimensional vector with data points (row-wise)
   * @param result the result vector of the matrix vector multiplication
   */
  void mult_transpose(sgpp::base::GridStorage& storage, BASIS& basis,
                      sgpp::base::DataVector& source, DataMatrix& x,
                      DataVectorDistributed& result) {
    if (!result.isProcessMapped()) {
      return;
    }

    result.setAll(0.0);
    size_t source_size = source.getSize();

    sgpp::base::DataVector privateResult(result.getGlobalRows());
    privateResult.setAll(0.0);

    sgpp::base::DataVector line(x.getNcols());
    sgpp::base::AlgorithmEvaluationTransposed<BASIS> AlgoEvalTrans(storage);

    privateResult.setAll(0.0);

    // distribute into approximately even blocks
    auto processGrid = result.getProcessGrid();
    int currentProcess = processGrid->getRowColumnIndex();
    int numberOfProcesses = processGrid->getProcessesInGrid();

    size_t blockSize = static_cast<size_t>(source_size / numberOfProcesses);

    size_t startIndex = currentProcess * blockSize;
    size_t endIndex = source_size;

    if (currentProcess != numberOfProcesses - 1) {
      endIndex = startIndex + blockSize;
    }

    for (size_t i = startIndex; i < endIndex; i++) {
      x.getRow(i, line);

      AlgoEvalTrans(basis, line, source[i], privateResult);
    }

    // distributed local results: for loop needed, loop through processes, distribute with different
    // master row/column needed
    DataVectorDistributed distributedResult(processGrid, result.getGlobalRows(),
                                            result.getBlockSize());
    for (int i = 0; i < processGrid->getTotalRows(); i++) {
      for (int j = 0; j < processGrid->getTotalColumns(); j++) {
        // distribute from this process as master and add partial result to total result
        distributedResult.distribute(privateResult.data(), i, j);
        result.add(distributedResult);
      }
    }
  }

  /**
   * Performs a distributed mass evaluation
   *
   * @param storage GridStorage object that contains the grid's points information
   * @param basis a reference to a class that implements a specific basis
   * @param source the coefficients of the grid points
   * @param x the d-dimensional vector with data points (row-wise)
   * @param result the result vector of the matrix vector multiplication
   */
  void mult(sgpp::base::GridStorage& storage, BASIS& basis, sgpp::base::DataVector& source,
            DataMatrix& x, DataVectorDistributed& result) {
    if (!result.isProcessMapped()) {
      return;
    }

    sgpp::base::DataVector line(x.getNcols());
    sgpp::base::AlgorithmEvaluation<BASIS> AlgoEval(storage);

    for (size_t i = 0; i < result.getLocalRows(); i++) {
      size_t globalRow = result.localToGlobalRowIndex(i);
      x.getRow(globalRow, line);

      result.getLocalPointer()[i] = AlgoEval(basis, line, source);
    }
  }
};

}  // namespace datadriven
}  // namespace sgpp
