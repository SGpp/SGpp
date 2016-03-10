// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALPREWAVELET_HPP
#define OPERATIONEVALPREWAVELET_HPP

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * This class implements OperationEval for a grids with prewavelet basis ansatzfunctions without boundaries
 */
class OperationEvalPrewavelet : public OperationEval {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   */
  explicit OperationEvalPrewavelet(GridStorage& storage) : storage(storage) {}

  /**
   * Destructor
   */
  ~OperationEvalPrewavelet() override {}

  double eval(const DataVector& alpha,
               const DataVector& point) override;
  virtual double test(const DataVector& alpha,
                       const DataVector& data,
                       const DataVector& classes);
  virtual double integrate(const DataVector& alpha);

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONEVELMODLINEAR_HPP */
