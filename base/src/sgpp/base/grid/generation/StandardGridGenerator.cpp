// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementInteraction.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace base {

StandardGridGenerator::StandardGridGenerator(GridStorage& storage) : storage(storage) {}

StandardGridGenerator::~StandardGridGenerator() {}

void StandardGridGenerator::regular(size_t level) {
  HashGenerator gen;
  gen.regular(this->storage, static_cast<level_t>(level));
}

void StandardGridGenerator::regular(size_t level, double T) {
  HashGenerator gen;
  gen.regular(this->storage, static_cast<level_t>(level), T);
}

void StandardGridGenerator::regularInter(size_t level,
                                         const std::vector<std::vector<size_t>>& terms, double T) {
  HashGenerator gen;
  gen.regularInter(this->storage, static_cast<level_t>(level), terms, T);
}

void StandardGridGenerator::cliques(size_t level, size_t clique_size) {
  HashGenerator gen;
  gen.cliques(this->storage, static_cast<level_t>(level), clique_size);
}

void StandardGridGenerator::cliques(size_t level, size_t clique_size, double T) {
  HashGenerator gen;
  gen.cliques(this->storage, static_cast<level_t>(level), clique_size, T);
}

void StandardGridGenerator::full(size_t level) {
  HashGenerator gen;
  gen.full(this->storage, static_cast<level_t>(level));
}

void StandardGridGenerator::anisotropicFull(std::vector<size_t> dimlevels) {
  HashGenerator gen;
  gen.anisotropicFull(this->storage, dimlevels);
}

void StandardGridGenerator::refine(RefinementFunctor& func, std::vector<size_t>* addedPoints) {
  HashRefinement refine;
  refine.free_refine(this->storage, func, addedPoints);
}

void StandardGridGenerator::refineInter(RefinementFunctor& func,
                                        const std::unordered_set<std::vector<bool>>& interactions) {
  auto refine = HashRefinementInteraction(interactions);
  refine.free_refine(this->storage, func);
}

void StandardGridGenerator::refineInter(RefinementFunctor& func,
                                        const std::vector<std::vector<size_t>>& interactions) {
  auto interset = std::unordered_set<std::vector<bool>>();
  const auto dim = storage.getDimension();
  for (const auto& interaction : interactions) {
    auto term = std::vector<bool>(dim, false);
    for (const auto i : interaction) {
      term[i] = true;
    }
    interset.insert(term);
  }
  refineInter(func, interset);
}

size_t StandardGridGenerator::getNumberOfRefinablePoints() {
  HashRefinement refine;
  return refine.getNumberOfRefinablePoints(this->storage);
}

void StandardGridGenerator::coarsen(CoarseningFunctor& func, DataVector& alpha) {
  HashCoarsening coarsen;
  coarsen.free_coarsen(this->storage, func, alpha);
}

void StandardGridGenerator::coarsenNFirstOnly(CoarseningFunctor& func, DataVector& alpha,
                                              size_t numFirstOnly) {
  HashCoarsening coarsen;
  coarsen.free_coarsen_NFirstOnly(this->storage, func, alpha, numFirstOnly);
}

size_t StandardGridGenerator::getNumberOfRemovablePoints() {
  HashCoarsening coarsen;
  return coarsen.getNumberOfRemovablePoints(this->storage);
}

void StandardGridGenerator::refineMaxLevel(RefinementFunctor& func, size_t maxLevel) {
  throw generation_exception("StandardGridGenerator::refineMaxLevel is not implemented");
}

size_t StandardGridGenerator::getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) {
  throw generation_exception(
      "StandardGridGenerator::getNumberOfRefinablePointsToMaxLevel "
      "is not implemented");
  return 0;
}

}  // namespace base
}  // namespace sgpp
