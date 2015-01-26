/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Florian Zipperle (florian.zipperle@tum.de)

#include <sgpp/base/grid/generation/PeriodicGridGenerator.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

  PeriodicGridGenerator::PeriodicGridGenerator(GridStorage* storage) : storage(storage) {
    }

  PeriodicGridGenerator::~PeriodicGridGenerator() {
    }

    void PeriodicGridGenerator::regular(int level) {
      if (level < 0) {
        throw generation_exception("Grid level value is negative");
      }

      HashGenerator gen;
      gen.regularWithPeriodicBoundaries(storage, level);
    }

    void PeriodicGridGenerator::full(int level) {
      if (level < 0) {
        throw generation_exception("Grid level value is negative");
      }

      throw generation_exception("PeriodicGridGenerator::full is not implemented");
    }

    void PeriodicGridGenerator::refine(RefinementFunctor* func) {
      throw generation_exception("PeriodicGridGenerator::refine is not implemented");
    }

    size_t PeriodicGridGenerator::getNumberOfRefinablePoints() {
      throw generation_exception("PeriodicGridGenerator::getNumberOfRefinablePoints is not implemented");
      return 0;
    }

    void PeriodicGridGenerator::cliques(int level, size_t clique_size) {
      if (level < 0) {
      throw generation_exception("Grid level value is negative");
      }

      HashGenerator gen;
      gen.cliques(this->storage, static_cast<HashGenerator::level_t>(level), clique_size);
    }

    void PeriodicGridGenerator::coarsen(CoarseningFunctor* func, DataVector* alpha) {
      throw generation_exception("PeriodicGridGenerator::coarsen is not implemented");
    }

    void PeriodicGridGenerator::coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly) {
      throw generation_exception("PeriodicGridGenerator::coarsenNFirstOnly is not implemented");
    }

    size_t PeriodicGridGenerator::getNumberOfRemovablePoints() {

      throw generation_exception("PeriodicGridGenerator::getNumberOfRemovablePoints is not implemented");
      return 0;
    }

    void PeriodicGridGenerator::refineMaxLevel(RefinementFunctor* func, int maxLevel) {
      if (maxLevel < 0) {
        throw generation_exception("Grid level value is negative");
      }

      throw generation_exception("PeriodicGridGenerator::refineMaxLevel is not implemented");
    }

    size_t PeriodicGridGenerator::getNumberOfRefinablePointsToMaxLevel(int maxLevel) {
      if (maxLevel < 0) {
        throw generation_exception("Grid level value is negative");
      }

      throw generation_exception("PeriodicGridGenerator::getNumberOfRefinablePointsToMaxLevel is not implemented");
      return 0;
    }

  }
}

