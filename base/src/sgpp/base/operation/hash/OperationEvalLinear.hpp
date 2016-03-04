// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALLINEAR_HPP
#define OPERATIONEVALLINEAR_HPP

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * This class implements OperationEval for a grids with linear basis ansatzfunctions without boundaries
 */
class OperationEvalLinear : public OperationEval {
 public:
  /**
   * Constructor of OperationEvalLinear
   *
   * @param storage Pointer to the grid's gridstorage obejct
   */
  explicit OperationEvalLinear(GridStorage& storage) : storage(storage) {}

  /**
   * Destructor
   */
  ~OperationEvalLinear() override {}

  double eval(const DataVector& alpha,
               const DataVector& point) override;

 protected:
  /// reference to the grid's GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONEVALLINEAR_HPP */
