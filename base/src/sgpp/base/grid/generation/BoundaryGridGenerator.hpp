// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TRUNCATEDBOUNDARYGRIDGENERATOR_HPP
#define TRUNCATEDBOUNDARYGRIDGENERATOR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
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
   * @param boundaryLevel 1 + how much levels the boundary is coarser than
   *                      the main axes, 0 means one level finer,
   *                      1 means same level,
   *                      2 means one level coarser, etc.
   */
  explicit BoundaryGridGenerator(GridStorage& storage, level_t boundaryLevel = 1);

  /**
   * Destructor
   */
  ~BoundaryGridGenerator() override;

  level_t getBoundaryLevel() const;
  void setBoundaryLevel(level_t boundaryLevel);

  void regular(size_t level) override;
  void regular(size_t level, double T) override;
  void cliques(size_t level, size_t clique_size) override;
  void full(size_t level) override;
  void refine(RefinementFunctor& func, std::vector<size_t>* addedPoints = 0) override;
  size_t getNumberOfRefinablePoints() override;

  void coarsen(CoarseningFunctor& func, DataVector& alpha) override;
  void coarsenNFirstOnly(CoarseningFunctor& func, DataVector& alpha, size_t numFirstOnly) override;
  size_t getNumberOfRemovablePoints() override;

  void refineMaxLevel(RefinementFunctor& func, size_t maxLevel) override;
  size_t getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) override;

 protected:
  /// reference to the grid's storage object
  GridStorage& storage;
  /// level at which the boundary points should be inserted
  level_t boundaryLevel;
};

}  // namespace base
}  // namespace sgpp

#endif /* TRUNCATEDBOUNDARYGRIDGENERATOR_HPP */
