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
#include "TruncatedBoundaryGridGenerator.hpp"


namespace SGPP {
  namespace base {

    TruncatedBoundaryGridGenerator::TruncatedBoundaryGridGenerator(GridStorage* storage) : storage(storage) {
    }

    TruncatedBoundaryGridGenerator::~TruncatedBoundaryGridGenerator() {
    }

    void TruncatedBoundaryGridGenerator::regular(int level) {
      if (level < 0) {
        throw generation_exception("Grid level value is negative");
      }

      HashGenerator gen;
      gen.regularWithBoundaries(this->storage, static_cast<HashGenerator::level_t>(level), true);
    }

    void TruncatedBoundaryGridGenerator::cliques(int level, size_t clique_size) {
      throw generation_exception("Method is not implemented");
    }

    void TruncatedBoundaryGridGenerator::full(int level) {
      if (level < 0) {
        throw generation_exception("Grid level value is negative");
      }

      HashGenerator gen;
      gen.fullWithBoundary(this->storage, static_cast<HashGenerator::level_t>(level));
    }

    void TruncatedBoundaryGridGenerator::refine(RefinementFunctor* func) {
      HashRefinementBoundaries refine;
      refine.free_refine(this->storage, func);
    }

    size_t TruncatedBoundaryGridGenerator::getNumberOfRefinablePoints() {
      HashRefinementBoundaries refine;
      return refine.getNumberOfRefinablePoints(this->storage);
    }

    void TruncatedBoundaryGridGenerator::coarsen(CoarseningFunctor* func, DataVector* alpha) {
      HashCoarsening coarsen;
      coarsen.free_coarsen(this->storage, func, alpha);
    }

    void TruncatedBoundaryGridGenerator::coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly) {
      HashCoarsening coarsen;
      coarsen.free_coarsen_NFirstOnly(this->storage, func, alpha, numFirstOnly);
    }

    size_t TruncatedBoundaryGridGenerator::getNumberOfRemovablePoints() {
      HashCoarsening coarsen;
      return coarsen.getNumberOfRemovablePoints(this->storage);
    }

    void TruncatedBoundaryGridGenerator::refineMaxLevel(RefinementFunctor* func, int maxLevel) {
      if (maxLevel < 0) {
        throw generation_exception("Grid level value is negative");
      }

      HashRefinementBoundariesMaxLevel refine;
      refine.refineToMaxLevel(this->storage, func, static_cast<HashGenerator::level_t>(maxLevel));
    }

    size_t TruncatedBoundaryGridGenerator::getNumberOfRefinablePointsToMaxLevel(int maxLevel) {
      if (maxLevel < 0) {
        throw generation_exception("Grid level value is negative");
      }

      HashRefinementBoundariesMaxLevel refine;
      return refine.getNumberOfRefinablePointsToMaxLevel(this->storage, static_cast<HashGenerator::level_t>(maxLevel));
    }

  }
}