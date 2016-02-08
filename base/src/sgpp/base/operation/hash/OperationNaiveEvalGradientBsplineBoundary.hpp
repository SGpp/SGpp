// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALGRADIENTBSPLINEBOUNDARY_HPP
#define OPERATIONEVALGRADIENTBSPLINEBOUNDARY_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradient.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBoundaryBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
namespace base {

/**
 * Operation for evaluating B-spline linear combinations on Boundary grids and their gradients.
 */
class OperationNaiveEvalGradientBsplineBoundary : public
  OperationNaiveEvalGradient {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   * @param degree    B-spline degree
   */
  OperationNaiveEvalGradientBsplineBoundary(GridStorage* storage, size_t degree) :
    storage(storage), base(degree) {
  }

  /**
   * Destructor.
   */
  virtual ~OperationNaiveEvalGradientBsplineBoundary() override {
  }

  /**
   * @param       alpha       coefficient vector
   * @param       point       evaluation point
   * @param[out]  gradient    gradient of linear combination
   * @return                  value of linear combination
   */
  virtual float_t evalGradient(const DataVector& alpha,
                               const DataVector& point,
                               DataVector& gradient) override;

 protected:
  /// storage of the sparse grid
  GridStorage* storage;
  /// 1D B-spline basis
  SBsplineBoundaryBase base;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONEVALGRADIENTBSPLINEBOUNDARY_HPP */
