// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVALBSPLINECLENSHAWCURTIS_HPP
#define OPERATIONNAIVEEVALBSPLINECLENSHAWCURTIS_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineClenshawCurtisBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
namespace base {

/**
 * Operation for evaluating B-spline linear combinations on Clenshaw-Curtis grids.
 */
class OperationNaiveEvalBsplineClenshawCurtis : public OperationNaiveEval {
 public:
  /**
   * Constructor.
   *
   * @param storage       storage of the sparse grid
   * @param degree        B-spline degree
   */
  OperationNaiveEvalBsplineClenshawCurtis(
    GridStorage& storage, size_t degree)
    : storage(storage),
      base(degree) {
  }

  /**
   * Destructor.
   */
  ~OperationNaiveEvalBsplineClenshawCurtis() override {
  }

  /**
   * @param alpha     coefficient vector
   * @param point     evaluation point
   * @return          value of linear combination
   */
  float_t eval(const DataVector& alpha, const DataVector& point) override;

 protected:
  /// storage of the sparse grid
  GridStorage& storage;
  /// 1D B-spline basis
  SBsplineClenshawCurtisBase base;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONNAIVEEVALBSPLINECLENSHAWCURTIS_HPP */
