// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundariesMaxLevel.hpp>
#include <sgpp/globaldef.hpp>
#include <vector>

namespace sgpp {
namespace base {

BoundaryGridGenerator::BoundaryGridGenerator(GridStorage& storage, level_t boundaryLevel)
    : storage(storage), boundaryLevel(boundaryLevel) {}

BoundaryGridGenerator::~BoundaryGridGenerator() {}

level_t BoundaryGridGenerator::getBoundaryLevel() const { return boundaryLevel; }

void BoundaryGridGenerator::setBoundaryLevel(level_t boundaryLevel) {
  this->boundaryLevel = boundaryLevel;
}

void BoundaryGridGenerator::regular(size_t level) {
  HashGenerator gen;
  gen.regularWithBoundaries(this->storage, static_cast<level_t>(level), boundaryLevel);
}

void BoundaryGridGenerator::cliques(size_t level, size_t clique_size) {
  throw generation_exception("Method is not implemented");
}

void BoundaryGridGenerator::full(size_t level) {
  HashGenerator gen;
  gen.fullWithBoundary(this->storage, static_cast<level_t>(level));
}

void BoundaryGridGenerator::refine(RefinementFunctor& func, std::vector<size_t>* addedPoints) {
  HashRefinementBoundaries refine;
  refine.free_refine(this->storage, func, addedPoints);
}

size_t BoundaryGridGenerator::getNumberOfRefinablePoints() {
  HashRefinementBoundaries refine;
  return refine.getNumberOfRefinablePoints(this->storage);
}

void BoundaryGridGenerator::coarsen(CoarseningFunctor& func, std::vector<size_t>* removedSeq) {
  HashCoarsening coarsen;
  coarsen.free_coarsen(this->storage, func, nullptr, removedSeq);
}

void BoundaryGridGenerator::coarsenNFirstOnly(CoarseningFunctor& func, size_t numFirstOnly,
                                              std::vector<size_t>* removedSeq,
                                              size_t minIndexConsidered) {
  HashCoarsening coarsen;
  coarsen.free_coarsen_NFirstOnly(this->storage, func, numFirstOnly, minIndexConsidered, nullptr,
                                  removedSeq);
}

size_t BoundaryGridGenerator::getNumberOfRemovablePoints() {
  HashCoarsening coarsen;
  return coarsen.getNumberOfRemovablePoints(this->storage);
}

void BoundaryGridGenerator::refineMaxLevel(RefinementFunctor& func, size_t maxLevel) {
  HashRefinementBoundariesMaxLevel refine;
  refine.refineToMaxLevel(this->storage, func, static_cast<level_t>(maxLevel));
}

size_t BoundaryGridGenerator::getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) {
  HashRefinementBoundariesMaxLevel refine;
  return refine.getNumberOfRefinablePointsToMaxLevel(this->storage, static_cast<level_t>(maxLevel));
}

}  // namespace base
}  // namespace sgpp
