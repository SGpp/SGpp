// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONQUADRATUREMODLINEAR_HPP
#define OPERATIONQUADRATUREMODLINEAR_HPP

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, mod linear grid
 */
class OperationQuadratureModLinear : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureModLinear
   *
   * @param storage Pointer to the grid's GridStorage object
   */
  explicit OperationQuadratureModLinear(GridStorage& storage) : storage(storage) {}

  ~OperationQuadratureModLinear() override {}
  double doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  // ModLinear Basis object
  SLinearModifiedBase base;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONQUADRATURE_HPP */
