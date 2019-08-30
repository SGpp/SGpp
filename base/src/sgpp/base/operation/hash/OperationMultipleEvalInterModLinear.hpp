// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALINTERMODLINEAR_HPP
#define OPERATIONMULTIPLEEVALINTERMODLINEAR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace base {

/**
 * This class implements OperationMultipleEval for a grids with mod linear basis ansatzfunctions
 * and a set of interactions that limit the subspaces that are included.
 *
 * The evaluation iterates over all allowed subspaces and only viewes ansatzfunctions of level > 1
 * since the level one ansatz function evaluates to 1
 *
 *
 * This Operation is more efficeint than the OperationMultipleEvalModLinear.
 */
class OperationMultipleEvalInterModLinear : public OperationMultipleEval {
 public:
  /**
   * Constructor
   *
   * @param grid grid
   * @param dataset the dataset that should be evaluated
   * @param interactions the interactionterms that limit what subspaces are included
   */
  OperationMultipleEvalInterModLinear(Grid& grid, DataMatrix& dataset,
    std::vector<std::vector<size_t>>& interactions)
      : OperationMultipleEval(grid, dataset), storage(grid.getStorage()) {
        this->interactions = interactions;
      }

  /**
   * Destructor
   */
  ~OperationMultipleEvalInterModLinear() override {}

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;


  double getDuration() override { return 0.0; }


 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
  // interactions
  std::vector<std::vector<size_t>> interactions;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONMULTIPLEEVALMODLINEAR_HPP */
