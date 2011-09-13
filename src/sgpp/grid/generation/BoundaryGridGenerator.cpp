/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/generation/BoundaryGridGenerator.hpp"
#include "grid/GridStorage.hpp"

#include "grid/generation/hashmap/HashRefinementBoundariesMaxLevel.hpp"
#include "grid/generation/hashmap/HashCoarsening.hpp"
#include "grid/generation/hashmap/HashRefinementBoundaries.hpp"
#include "grid/generation/hashmap/HashGenerator.hpp"

namespace sg
{
namespace base
{

BoundaryGridGenerator::BoundaryGridGenerator(GridStorage* storage) : storage(storage)
{
}

BoundaryGridGenerator::~BoundaryGridGenerator()
{
}

void BoundaryGridGenerator::regular(size_t level)
{
	HashGenerator gen;
	gen.regularWithBoundaries(this->storage, level, false);
}

void BoundaryGridGenerator::refine(RefinementFunctor* func)
{
	HashRefinementBoundaries refine;
	refine.free_refine(this->storage, func);
}

size_t BoundaryGridGenerator::getNumberOfRefinablePoints()
{
	HashRefinementBoundaries refine;
	return refine.getNumberOfRefinablePoints(this->storage);
}

void BoundaryGridGenerator::coarsen(CoarseningFunctor* func, DataVector* alpha)
{
	HashCoarsening coarsen;
	coarsen.free_coarsen(this->storage, func, alpha);
}

void BoundaryGridGenerator::coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly)
{
	HashCoarsening coarsen;
	coarsen.free_coarsen_NFirstOnly(this->storage, func, alpha, numFirstOnly);
}

size_t BoundaryGridGenerator::getNumberOfRemoveablePoints()
{
	HashCoarsening coarsen;
	return coarsen.getNumberOfRemovablePoints(this->storage);
}

void BoundaryGridGenerator::refineMaxLevel(RefinementFunctor* func, unsigned int maxLevel)
{
	HashRefinementBoundariesMaxLevel refine;
	refine.refineToMaxLevel(this->storage, func, maxLevel);
}

size_t BoundaryGridGenerator::getNumberOfRefinablePointsToMaxLevel(unsigned int maxLevel)
{
	HashRefinementBoundariesMaxLevel refine;
	return refine.getNumberOfRefinablePointsToMaxLevel(this->storage, maxLevel);
}

}
}
