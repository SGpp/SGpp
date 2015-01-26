/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/

#include <sgpp/base/grid/generation/TruncatedTrapezoidGridGenerator.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    TruncatedTrapezoidGridGenerator::TruncatedTrapezoidGridGenerator(GridStorage* storage) : storage(storage) {
    }

    TruncatedTrapezoidGridGenerator::~TruncatedTrapezoidGridGenerator() {
    }

    void TruncatedTrapezoidGridGenerator::regular(int level) {
      if (level < 0) {
        throw generation_exception("Grid level value is negative");
      }

      this->truncated( static_cast<HashGenerator::level_t>(level), 1);
    }

    void TruncatedTrapezoidGridGenerator::cliques(int level, size_t clique_size) {
		throw generation_exception("Method is not implemented");
	}


    void TruncatedTrapezoidGridGenerator::truncated(int level, int l_user) {
      if (level < 0) {
        throw generation_exception("Grid level value is negative");
      }

      HashGenerator gen;
      gen.truncated(this->storage, static_cast<HashGenerator::level_t>(level), static_cast<HashGenerator::level_t>(l_user));
    }

    //void BoundaryGridGenerator::refine(RefinementFunctor* func)
    //{
    //  HashRefinementBoundaries refine;
    //  refine.free_refine(this->storage, func);
    //}
    //
    //size_t BoundaryGridGenerator::getNumberOfRefinablePoints()
    //{
    //  HashRefinementBoundaries refine;
    //  return refine.getNumberOfRefinablePoints(this->storage);
    //}
    //
    //void BoundaryGridGenerator::coarsen(CoarseningFunctor* func, DataVector* alpha)
    //{
    //  HashCoarsening coarsen;
    //  coarsen.free_coarsen(this->storage, func, alpha);
    //}
    //
    //size_t BoundaryGridGenerator::getNumberOfRemovablePoints()
    //{
    //  HashCoarsening coarsen;
    //  return coarsen.getNumberOfRemovablePoints(this->storage);
    //}
    //
    //void BoundaryGridGenerator::refineMaxLevel(RefinementFunctor* func, unsigned int maxLevel)
    //{
    //  HashRefinementBoundariesMaxLevel refine;
    //  refine.refineToMaxLevel(this->storage, func, maxLevel);
    //}
    //
    //size_t BoundaryGridGenerator::getNumberOfRefinablePointsToMaxLevel(unsigned int maxLevel)
    //{
    //  HashRefinementBoundariesMaxLevel refine;
    //  return refine.getNumberOfRefinablePointsToMaxLevel(this->storage, maxLevel);
    //}

  }
}


