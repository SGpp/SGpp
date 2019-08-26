// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONSTENCILHIERARCHISATIONLINEAR_HPP
#define OPERATIONSTENCILHIERARCHISATIONLINEAR_HPP

#include <sgpp/base/operation/hash/OperationStencilHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>


namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, linear grid without boundaries
 */
class OperationStencilHierarchisationLinear : public
  OperationStencilHierarchisation {
 public:
  /**
   * Constructor of OperationStencilHierarchisationLinear
   *
   * @param storage Pointer to the grid's gridstorage obejct
   */
  explicit OperationStencilHierarchisationLinear(GridStorage& storage) :
    storage(storage),
    surplusStencil(0), neighborStencil(0), weightStencil(0) {}

  /**
   * Destructor
   */
  ~OperationStencilHierarchisationLinear() override {}

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;


  const IndexStencil& getSurplusStencil() const override {
    return surplusStencil;
  }

  const IndexStencil& getNeighborStencil() const override {
    return neighborStencil;
  }

  const WeightStencil& getWeightStencil() const override {
    return weightStencil;
  }

  size_t getStencilSize() const override {
    return surplusStencil.size();
  }

 protected:
  /// reference to the grid's GridStorage object
  GridStorage& storage;

  /// Index array with surplus indices
  IndexStencil surplusStencil;

  /// Index array with neighboring surplus indices
  IndexStencil neighborStencil;

  /// Index array with surplus indices
  WeightStencil weightStencil;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONSTENCILHIERARCHISATION_HPP */
