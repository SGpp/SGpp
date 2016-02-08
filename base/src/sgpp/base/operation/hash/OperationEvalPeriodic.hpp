// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALPERIODIC_HPP
#define OPERATIONEVALPERIODIC_HPP

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * This class implements OperationEval for a grids with periodic linear basis ansatzfunctions with
 *
 */
class OperationEvalPeriodic : public OperationEval {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   */
  OperationEvalPeriodic(GridStorage* storage) : storage(storage) {}

  /**
   * Destructor
   */
  virtual ~OperationEvalPeriodic() override {}

  virtual float_t eval(const DataVector& alpha,
                       const DataVector& point) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage* storage;

};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONEVALPERIODIC_HPP */
