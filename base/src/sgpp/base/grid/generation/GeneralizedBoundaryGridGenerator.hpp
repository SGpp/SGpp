// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GENERALIZEDTRUNCATEDBOUNDARYGRIDGENERATOR_HPP_
#define GENERALIZEDTRUNCATEDBOUNDARYGRIDGENERATOR_HPP_

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
class GeneralizedBoundaryGridGenerator : public GridGenerator {
 public:
  /**
   * Constructor
   *
   * @param storage template type that holds the grid points
   */
  explicit GeneralizedBoundaryGridGenerator(GridStorage& storage);

  /**
   * Destructor
   */
  ~GeneralizedBoundaryGridGenerator() override;
  /**
   * Creates a regular truncated boundary grid with given level and l_user=1
   * Is the same as the regular truncated grid
   * */
  void regular(size_t level) override;
  void cliques(size_t level, size_t clique_size) override;
  void full(size_t level) override {}
  /**
   * Creates a super truncated boundary grid with given level and l_user
   * @param level the maximum level of the grid
   * @param l_user the number of fullgrids cut off from the boundaries.
   * */
  void truncated(size_t level, size_t l_user) override;
  void refine(RefinementFunctor& func, std::vector<size_t>* addedPoints = nullptr) override {}
  size_t getNumberOfRefinablePoints() override { return 0; }

  void coarsen(CoarseningFunctor& func, DataVector& alpha,
               std::vector<size_t>* removedSeq) override {}
  void coarsenNFirstOnly(CoarseningFunctor& func, DataVector& alpha, size_t numFirstOnly,
                         std::vector<size_t>* removedSeq) override {}
  size_t getNumberOfRemovablePoints() override { return 0; }

  void refineMaxLevel(RefinementFunctor& func, size_t maxLevel) override {
    throw generation_exception("refineMaxLevel is not implemented");
  }
  size_t getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) override {
    throw generation_exception("getNumberOfRefinablePointsToMaxLevel is not implemented");
    return 0;
  }

 protected:
  /// reference to the grid's storage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif /* GENERALIZEDTRUNCATEDBOUNDARYGRIDGENERATOR_HPP_ */
