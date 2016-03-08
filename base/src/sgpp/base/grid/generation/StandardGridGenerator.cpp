// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

StandardGridGenerator::StandardGridGenerator(GridStorage& storage) : storage(
    storage) {
}

StandardGridGenerator::~StandardGridGenerator() {
}

void StandardGridGenerator::regular(size_t level) {
  HashGenerator gen;
  gen.regular(this->storage, static_cast<HashGenerator::level_t>(level));
}

void StandardGridGenerator::cliques(size_t level, size_t clique_size) {
  HashGenerator gen;
  gen.cliques(this->storage, static_cast<HashGenerator::level_t>(level),
              clique_size);
}

void StandardGridGenerator::full(size_t level) {
  HashGenerator gen;
  gen.full(this->storage, static_cast<HashGenerator::level_t>(level));
}

void StandardGridGenerator::refine(RefinementFunctor& func) {
  HashRefinement refine;
  refine.free_refine(this->storage, func);
}

size_t StandardGridGenerator::getNumberOfRefinablePoints() {
  HashRefinement refine;
  return refine.getNumberOfRefinablePoints(this->storage);
}

void StandardGridGenerator::coarsen(CoarseningFunctor& func,
                                    DataVector& alpha) {
  HashCoarsening coarsen;
  coarsen.free_coarsen(this->storage, func, alpha);
}

void StandardGridGenerator::coarsenNFirstOnly(CoarseningFunctor& func,
    DataVector& alpha, size_t numFirstOnly) {
  HashCoarsening coarsen;
  coarsen.free_coarsen_NFirstOnly(this->storage, func, alpha, numFirstOnly);
}

size_t StandardGridGenerator::getNumberOfRemovablePoints() {
  HashCoarsening coarsen;
  return coarsen.getNumberOfRemovablePoints(this->storage);
}

void StandardGridGenerator::refineMaxLevel(RefinementFunctor& func,
    size_t maxLevel) {
  throw generation_exception(
    "StandardGridGenerator::refineMaxLevel is not implemented");
}

size_t StandardGridGenerator::getNumberOfRefinablePointsToMaxLevel(
  size_t maxLevel) {
  throw generation_exception(
    "StandardGridGenerator::getNumberOfRefinablePointsToMaxLevel "
    "is not implemented");
  return 0;
}

}  // namespace base
}  // namespace sgpp
