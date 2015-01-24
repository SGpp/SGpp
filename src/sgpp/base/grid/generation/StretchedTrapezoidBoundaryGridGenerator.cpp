/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "base/grid/generation/StretchedTrapezoidBoundaryGridGenerator.hpp"
#include "base/grid/GridStorage.hpp"

#include "base/grid/generation/hashmap/HashRefinementBoundariesMaxLevel.hpp"
#include "base/grid/generation/hashmap/HashCoarsening.hpp"
#include "base/grid/generation/hashmap/HashRefinementBoundaries.hpp"
#include "base/grid/generation/hashmap/HashGenerator.hpp"

namespace sg {
  namespace base {

    StretchedTrapezoidBoundaryGridGenerator::StretchedTrapezoidBoundaryGridGenerator(GridStorage* storage) : storage(storage) {
    }

    StretchedTrapezoidBoundaryGridGenerator::~StretchedTrapezoidBoundaryGridGenerator() {
    }

    void StretchedTrapezoidBoundaryGridGenerator::regular(int level) {
      if (level < 0) {
        throw generation_exception("Grid level value is negative");
      }

      HashGenerator gen;
      gen.regularWithBoundaries(this->storage, static_cast<HashGenerator::level_t>(level), true);
    }

    void StretchedTrapezoidBoundaryGridGenerator::cliques(int level, size_t clique_size) {
		throw generation_exception("Method is not implemented");
	}

    void StretchedTrapezoidBoundaryGridGenerator::full(int level) {
      if (level < 0) {
        throw generation_exception("Grid level value is negative");
      }

      HashGenerator gen;
      gen.fullWithBoundary(this->storage, static_cast<HashGenerator::level_t>(level));
    }

    void StretchedTrapezoidBoundaryGridGenerator::refine(RefinementFunctor* func) {
      HashRefinementBoundaries refine;
      refine.free_refine(this->storage, func);
    }

    size_t StretchedTrapezoidBoundaryGridGenerator::getNumberOfRefinablePoints() {
      HashRefinementBoundaries refine;
      return refine.getNumberOfRefinablePoints(this->storage);
    }

    void StretchedTrapezoidBoundaryGridGenerator::coarsen(CoarseningFunctor* func, DataVector* alpha) {
      HashCoarsening coarsen;
      coarsen.free_coarsen(this->storage, func, alpha);
    }

    void StretchedTrapezoidBoundaryGridGenerator::coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly) {
      HashCoarsening coarsen;
      coarsen.free_coarsen_NFirstOnly(this->storage, func, alpha, numFirstOnly);
    }

    size_t StretchedTrapezoidBoundaryGridGenerator::getNumberOfRemovablePoints() {
      HashCoarsening coarsen;
      return coarsen.getNumberOfRemovablePoints(this->storage);
    }

    void StretchedTrapezoidBoundaryGridGenerator::refineMaxLevel(RefinementFunctor* func, int maxLevel) {
      if (maxLevel < 0) {
        throw generation_exception("Grid level value is negative");
      }

      HashRefinementBoundariesMaxLevel refine;
      refine.refineToMaxLevel(this->storage, func, static_cast<HashGenerator::level_t>(maxLevel));
    }

    size_t StretchedTrapezoidBoundaryGridGenerator::getNumberOfRefinablePointsToMaxLevel(int maxLevel) {
      if (maxLevel < 0) {
        throw generation_exception("Grid level value is negative");
      }

      HashRefinementBoundariesMaxLevel refine;
      return refine.getNumberOfRefinablePointsToMaxLevel(this->storage, static_cast<HashGenerator::level_t>(maxLevel));
    }

  }
}
