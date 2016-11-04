// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * This class implements OperationEval for a grids with linear basis ansatzfunctions without
 * boundaries
 */
class OperationEvalLinearClenshawCurtis : public OperationEval {
 public:
  /**
   * Constructor of OperationEvalLinearClenshawCurtis
   *
   * @param storage Pointer to the grid's gridstorage obejct
   */
  explicit OperationEvalLinearClenshawCurtis(GridStorage& storage) : storage(storage) {}

  /**
   * Destructor
   */
  ~OperationEvalLinearClenshawCurtis() override {}

  double eval(const DataVector& alpha, const DataVector& point) override;

 protected:
  /// reference to the grid's GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp
