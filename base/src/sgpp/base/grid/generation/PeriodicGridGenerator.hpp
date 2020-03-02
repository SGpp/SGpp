// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PERIODICGRIDGENERATOR_HPP
#define PERIODICGRIDGENERATOR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/globaldef.hpp>
#include <vector>

namespace sgpp {
namespace base {

/**
 * GridGenerator for periodic grids with boundaries
 */
class PeriodicGridGenerator : public GridGenerator {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's storage object
   */
  explicit PeriodicGridGenerator(GridStorage& storage);

  /**
   * Destructor
   */
  ~PeriodicGridGenerator() override;

  void regular(size_t level) override;
  void regular(size_t level, double T) override;
  void full(size_t level) override;
  void refine(RefinementFunctor& func, std::vector<size_t>* addedPoints = nullptr) override;
  void cliques(size_t level, size_t clique_size) override;
  void cliques(size_t level, size_t clique_size, double T) override;
  size_t getNumberOfRefinablePoints() override;

  void coarsen(CoarseningFunctor& func, std::vector<size_t>* removedSeq) override;
  void coarsenNFirstOnly(CoarseningFunctor& func, size_t numFirstOnly,
                         std::vector<size_t>* removedSeq, size_t minIndexConsidered) override;
  size_t getNumberOfRemovablePoints() override;

  void refineMaxLevel(RefinementFunctor& func, size_t maxLevel) override;
  size_t getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) override;

 protected:
  /// reference to the storage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif /* PERIODICGRIDGENERATOR_HPP */
