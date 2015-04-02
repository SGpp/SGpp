// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundariesMaxLevel.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>

#include <sgpp/globaldef.hpp>
#include "StretchedTruncatedBoundaryGridGenerator.hpp"


namespace SGPP {
  namespace base {

    StretchedTruncatedBoundaryGridGenerator::StretchedTruncatedBoundaryGridGenerator(GridStorage* storage) : storage(storage) {
    }

    StretchedTruncatedBoundaryGridGenerator::~StretchedTruncatedBoundaryGridGenerator() {
    }

    void StretchedTruncatedBoundaryGridGenerator::regular(size_t level) {
      HashGenerator gen;
      gen.regularWithBoundaries(this->storage, static_cast<HashGenerator::level_t>(level), true);
    }

    void StretchedTruncatedBoundaryGridGenerator::cliques(size_t level, size_t clique_size) {
      throw generation_exception("Method is not implemented");
    }

    void StretchedTruncatedBoundaryGridGenerator::full(size_t level) {
      HashGenerator gen;
      gen.fullWithBoundary(this->storage, static_cast<HashGenerator::level_t>(level));
    }

    void StretchedTruncatedBoundaryGridGenerator::refine(RefinementFunctor* func) {
      HashRefinementBoundaries refine;
      refine.free_refine(this->storage, func);
    }

    size_t StretchedTruncatedBoundaryGridGenerator::getNumberOfRefinablePoints() {
      HashRefinementBoundaries refine;
      return refine.getNumberOfRefinablePoints(this->storage);
    }

    void StretchedTruncatedBoundaryGridGenerator::coarsen(CoarseningFunctor* func, DataVector* alpha) {
      HashCoarsening coarsen;
      coarsen.free_coarsen(this->storage, func, alpha);
    }

    void StretchedTruncatedBoundaryGridGenerator::coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly) {
      HashCoarsening coarsen;
      coarsen.free_coarsen_NFirstOnly(this->storage, func, alpha, numFirstOnly);
    }

    size_t StretchedTruncatedBoundaryGridGenerator::getNumberOfRemovablePoints() {
      HashCoarsening coarsen;
      return coarsen.getNumberOfRemovablePoints(this->storage);
    }

    void StretchedTruncatedBoundaryGridGenerator::refineMaxLevel(RefinementFunctor* func, size_t maxLevel) {
      HashRefinementBoundariesMaxLevel refine;
      refine.refineToMaxLevel(this->storage, func, static_cast<HashGenerator::level_t>(maxLevel));
    }

    size_t StretchedTruncatedBoundaryGridGenerator::getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) {
      HashRefinementBoundariesMaxLevel refine;
      return refine.getNumberOfRefinablePointsToMaxLevel(this->storage, static_cast<HashGenerator::level_t>(maxLevel));
    }

  }
}
