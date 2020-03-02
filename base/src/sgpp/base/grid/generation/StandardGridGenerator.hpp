// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STANDARDGRIDGENERATOR_HPP
#define STANDARDGRIDGENERATOR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/globaldef.hpp>

#include <set>
#include <unordered_set>
#include <vector>

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
  void regular(size_t level, double T) override;
  void regularInter(size_t level, const std::set<std::set<size_t>>& terms,
                    double T) override;
  void cliques(size_t level, size_t clique_size) override;
  void cliques(size_t level, size_t clique_size, double T) override;
  void full(size_t level) override;
  void anisotropicFull(std::vector<size_t> dimlevels) override;
  void refine(RefinementFunctor& func,
              std::vector<size_t>* addedPoints = nullptr) override;
  void refineInter(RefinementFunctor& func,
                   const std::unordered_set<std::vector<bool>>& interactions);
  void refineInter(RefinementFunctor& func,
                   const std::set<std::set<size_t>>& interactions) override;
  size_t getNumberOfRefinablePoints() override;

  void coarsen(CoarseningFunctor& func,
               std::vector<size_t>* removedSeq) override;
  void coarsenNFirstOnly(CoarseningFunctor& func, size_t numFirstOnly,
                         std::vector<size_t>* removedSeq,
                         size_t minIndexConsidered) override;
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
