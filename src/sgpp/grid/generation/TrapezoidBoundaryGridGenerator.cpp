/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/generation/TrapezoidBoundaryGridGenerator.hpp"
#include "grid/GridStorage.hpp"

#include "sgpp.hpp"

namespace sg
{

TrapezoidBoundaryGridGenerator::TrapezoidBoundaryGridGenerator(GridStorage* storage) : storage(storage)
{
}

TrapezoidBoundaryGridGenerator::~TrapezoidBoundaryGridGenerator()
{
}

void TrapezoidBoundaryGridGenerator::regular(size_t level)
{
	HashGenerator gen;
	gen.regularWithBoundaries(this->storage, level, true);
}

void TrapezoidBoundaryGridGenerator::refine(RefinementFunctor* func)
{
	HashRefinementBoundaries refine;
	refine.free_refine(this->storage, func);
}

size_t TrapezoidBoundaryGridGenerator::getNumberOfRefinablePoints()
{
	HashRefinementBoundaries refine;
	return refine.getNumberOfRefinablePoints(this->storage);
}

void TrapezoidBoundaryGridGenerator::coarsen(CoarseningFunctor* func, DataVector* alpha)
{
	HashCoarsening coarsen;
	coarsen.free_coarsen(this->storage, func, alpha);
}

size_t TrapezoidBoundaryGridGenerator::getNumberOfRemoveablePoints()
{
	HashCoarsening coarsen;
	return coarsen.getNumberOfRemovablePoints(this->storage);
}

void TrapezoidBoundaryGridGenerator::refineMaxLevel(RefinementFunctor* func, unsigned int maxLevel)
{
	HashRefinementBoundariesMaxLevel refine;
	refine.refineToMaxLevel(this->storage, func, maxLevel);
}

size_t TrapezoidBoundaryGridGenerator::getNumberOfRefinablePointsToMaxLevel(unsigned int maxLevel)
{
	HashRefinementBoundariesMaxLevel refine;
	return refine.getNumberOfRefinablePointsToMaxLevel(this->storage, maxLevel);
}

}
