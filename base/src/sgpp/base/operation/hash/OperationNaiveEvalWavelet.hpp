// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVALWAVELET_HPP
#define OPERATIONNAIVEEVALWAVELET_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
namespace base {

/**
 * Operation for evaluating wavelet linear combinations on Noboundary grids.
 */
class OperationNaiveEvalWavelet : public OperationNaiveEval {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   */
  explicit OperationNaiveEvalWavelet(GridStorage* storage) : storage(storage) {
  }

  /**
   * Destructor.
   */
  ~OperationNaiveEvalWavelet() override {
  }

  /**
   * @param alpha     coefficient vector
   * @param point     evaluation point
   * @return          value of linear combination
   */
  float_t eval(const DataVector& alpha, const DataVector& point) override;

 protected:
  /// storage of the sparse grid
  GridStorage* storage;
  /// 1D wavelet basis
  SWaveletBase base;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONNAIVEEVALWAVELET_HPP */
