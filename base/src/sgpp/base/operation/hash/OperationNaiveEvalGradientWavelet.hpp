// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVALGRADIENTWAVELET_HPP
#define OPERATIONNAIVEEVALGRADIENTWAVELET_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradient.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
namespace base {

/**
 * Operation for evaluating wavelet linear combinations on Noboundary grids and their gradients.
 */
class OperationNaiveEvalGradientWavelet : public OperationNaiveEvalGradient {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   */
  explicit OperationNaiveEvalGradientWavelet(GridStorage* storage) : storage(
      storage) {
  }

  /**
   * Destructor.
   */
  ~OperationNaiveEvalGradientWavelet() override {
  }

  /**
   * @param       alpha       coefficient vector
   * @param       point       evaluation point
   * @param[out]  gradient    gradient of linear combination
   * @return                  value of linear combination
   */
  float_t evalGradient(const DataVector& alpha,
                       const DataVector& point,
                       DataVector& gradient) override;

 protected:
  /// storage of the sparse grid
  GridStorage* storage;
  /// 1D wavelet basis
  SWaveletBase base;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONNAIVEEVALGRADIENTWAVELET_HPP */
