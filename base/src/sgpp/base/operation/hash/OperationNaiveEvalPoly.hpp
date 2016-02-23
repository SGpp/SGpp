// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVALPOLY_HPP_
#define OPERATIONNAIVEEVALPOLY_HPP_

#include <sgpp/base/operation/hash/OperationNaiveEval.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
namespace base {

class OperationNaiveEvalPoly: public OperationNaiveEval {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   * @param degree    polynomial degree
   */
  OperationNaiveEvalPoly(GridStorage& storage, size_t degree) :
    storage(storage), base(degree) {
  }

  ~OperationNaiveEvalPoly() override {
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
  SPolyBase base;
};

}  // namespace base
}  // namespace SGPP

#endif /* BASE_SRC_SGPP_BASE_OPERATION_HASH_OPERATIONNAIVEEVALPOLY_HPP_ */
