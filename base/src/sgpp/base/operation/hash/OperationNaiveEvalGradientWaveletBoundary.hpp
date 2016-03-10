// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVALGRADIENTWAVELETBOUNDARY_HPP
#define OPERATIONNAIVEEVALGRADIENTWAVELETBOUNDARY_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradient.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace base {

/**
 * Operation for evaluating wavelet linear combinations on Boundary grids and their gradients.
 */
class OperationNaiveEvalGradientWaveletBoundary : public
  OperationNaiveEvalGradient {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   */
  explicit OperationNaiveEvalGradientWaveletBoundary(GridStorage& storage) :
    storage(
      storage) {
  }

  /**
   * Destructor.
   */
  ~OperationNaiveEvalGradientWaveletBoundary() override {
  }

  /**
   * @param       alpha       coefficient vector
   * @param       point       evaluation point
   * @param[out]  gradient    gradient of linear combination
   * @return                  value of linear combination
   */
  double evalGradient(const DataVector& alpha,
                       const DataVector& point,
                       DataVector& gradient) override;

 protected:
  /// storage of the sparse grid
  GridStorage& storage;
  /// 1D wavelet basis
  SWaveletBoundaryBase base;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONNAIVEEVALGRADIENTWAVELETBOUNDARY_HPP */
