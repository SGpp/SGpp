// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STANDARDGRIDGENERATOR_HPP
#define STANDARDGRIDGENERATOR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * GridGenerator for standard grids without boundaries
 */
class StandardGridGenerator : public GridGenerator {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's storage object
   */
  explicit StandardGridGenerator(GridStorage& storage);

  /**
   * Destructor
   */
  ~StandardGridGenerator() override;

  void regular(size_t level) override;
  void cliques(size_t level, size_t clique_size) override;
  void full(size_t level) override;
  void refine(RefinementFunctor& func) override;
  size_t getNumberOfRefinablePoints() override;

  void coarsen(CoarseningFunctor& func, DataVector& alpha) override;
  void coarsenNFirstOnly(CoarseningFunctor& func, DataVector& alpha,
                         size_t numFirstOnly) override;
  size_t getNumberOfRemovablePoints() override;

  void refineMaxLevel(RefinementFunctor& func, size_t maxLevel) override;
  size_t getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) override;

 protected:
  /// reference to the storage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif /* STANDARDGRIDGEMERATOR_HPP */
