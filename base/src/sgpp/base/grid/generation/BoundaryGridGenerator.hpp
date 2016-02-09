// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TRUNCATEDBOUNDARYGRIDGENERATOR_HPP
#define TRUNCATEDBOUNDARYGRIDGENERATOR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * This class provides the interface for the grid generation
 * for grids with boundaries, pentagon cut through sub space scheme
 */
class BoundaryGridGenerator : public GridGenerator {
 public:
  /**
   * Constructor
   *
   * @param storage       template type that holds the grid points
   * @param boundaryLevel level at which the boundary points should be
   *                      inserted (default = 1: boundary has same level
   *                      as main axes)
   */
  explicit BoundaryGridGenerator(GridStorage* storage,
                                 level_t boundaryLevel = 1);

  /**
   * Destructor
   */
  ~BoundaryGridGenerator() override;

  void regular(size_t level) override;
  void cliques(size_t level, size_t clique_size) override;
  void full(size_t level) override;
  void refine(RefinementFunctor* func) override;
  size_t getNumberOfRefinablePoints() override;

  void coarsen(CoarseningFunctor* func, DataVector* alpha) override;
  void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha,
                         size_t numFirstOnly) override;
  size_t getNumberOfRemovablePoints() override;

  void refineMaxLevel(RefinementFunctor* func, size_t maxLevel) override;
  size_t getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) override;

 protected:
  /// Pointer to the grid's storage object
  GridStorage* storage;
  /// level at which the boundary points should be inserted
  level_t boundaryLevel;
};

}  // namespace base
}  // namespace SGPP

#endif /* TRUNCATEDBOUNDARYGRIDGENERATOR_HPP */
