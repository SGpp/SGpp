// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ALGORITHMMULTIPLEEVALUATION_HPP
#define ALGORITHMMULTIPLEEVALUATION_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/algorithm/AlgorithmEvaluation.hpp>
#include <sgpp/base/algorithm/AlgorithmEvaluationTransposed.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>
#include <utility>

namespace sgpp {
namespace base {

/**
 * Abstract implementation for multiple function evaluations. In Data Mining
 * to operators are needed: mass evaluation and transposed evaluation, referenced
 * in literature as matrix vector products with matrices B^T (mass evaluation) and B
 * (transposed evaluation).
 *
 * If there are @f$N@f$ basis functions @f$\varphi(\vec{x})@f$ and @f$m@f$ data points, then B is a
 * (mxN) matrix, with
 * @f[ (B)_{j,i} = \varphi_i(x_j). @f]
 *
 */
template <class BASIS>
class AlgorithmMultipleEvaluation {
 public:
  /**
   * Performs a transposed mass evaluation
   *
   * @param storage GridStorage object that contains the grid's points information
   * @param basis a reference to a class that implements a specific basis
   * @param source the coefficients of the grid points
   * @param x the d-dimensional vector with data points (row-wise)
   * @param result the result vector of the matrix vector multiplication
   */
  void mult_transpose(GridStorage& storage, BASIS& basis, DataVector& source, DataMatrix& x,
                      DataVector& result) {
    result.setAll(0.0);
    size_t source_size = source.getSize();

#pragma omp parallel
    {
      DataVector privateResult(result.getSize());
      privateResult.setAll(0.0);

      DataVector line(x.getNcols());
      AlgorithmEvaluationTransposed<BASIS> AlgoEvalTrans(storage);

      privateResult.setAll(0.0);

#pragma omp for schedule(static)

      for (size_t i = 0; i < source_size; i++) {
        x.getRow(i, line);

        AlgoEvalTrans(basis, line, source[i], privateResult);
      }

#pragma omp critical
      { result.add(privateResult); }
    }
  }
  // implementation requires OpenMP 4.0 support
  //        void mult_transpose(GridStorage& storage, BASIS& basis, DataVector& source,
  //                            DataMatrix& x, DataVector& result) {
  //          result.setAll(0.0);
  //          size_t source_size = source.getSize();
  //
  //          #pragma omp declare reduction(accumulate :
  // ...      sgpp::base::DataVector : omp_out.add(omp_in))
  // ...      initializer ( omp_priv = DataVector(omp_orig.getSize(), 0))
  //
  //          #pragma omp parallel
  //          {
  //            DataVector line(x.getNcols());
  //            AlgorithmEvaluationTransposed<BASIS> AlgoEvalTrans(storage);
  //
  //            #pragma omp for reduction(accumulate:result) schedule(static)
  //
  //            for (size_t i = 0; i < source_size; i++) {
  //              x.getRow(i, line);
  //
  //              AlgoEvalTrans(basis, line, source[i], result);
  //            }
  //          }
  //        }

  /**
   * Performs a mass evaluation
   *
   * @param storage GridStorage object that contains the grid's points information
   * @param basis a reference to a class that implements a specific basis
   * @param source the coefficients of the grid points
   * @param x the d-dimensional vector with data points (row-wise)
   * @param result the result vector of the matrix vector multiplication
   */
  void mult(GridStorage& storage, BASIS& basis, DataVector& source, DataMatrix& x,
            DataVector& result) {
    result.setAll(0.0);
    size_t result_size = result.getSize();

#pragma omp parallel
    {
      DataVector line(x.getNcols());
      AlgorithmEvaluation<BASIS> AlgoEval(storage);

#pragma omp for schedule(static)

      for (size_t i = 0; i < result_size; i++) {
        x.getRow(i, line);

        result[i] = AlgoEval(basis, line, source);
      }
    }
  }
};

}  // namespace base
}  // namespace sgpp

#endif /* ALGORITHMMULTIPLEEVALUATION_HPP */
