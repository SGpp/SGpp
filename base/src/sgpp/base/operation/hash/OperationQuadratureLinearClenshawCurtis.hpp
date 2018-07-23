// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, linear grid without boundaries
 */
class OperationQuadratureLinearClenshawCurtis : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureLinearClenshawCurtis
   *
   * @param storage Pointer to the grid's GridStorage object
   */
  explicit OperationQuadratureLinearClenshawCurtis(GridStorage& storage) : storage(storage) {}

  ~OperationQuadratureLinearClenshawCurtis() override {}

  /**
   * Quadrature for piecewise linear hat basis functions. Computes
   * @f[ \sum_{\vec{l}} 2^{-|\vec{l}|}\alpha_{\vec{l}}. @f]
   *
   * @param alpha Coefficient vector for current grid
   */
  double doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  /// Poly Boundary Basis object
  SLinearClenshawCurtisBase base;
};

}  // namespace base
}  // namespace sgpp
