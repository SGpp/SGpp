/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include "grid/generation/StretchedTrapezoidBoundaryGridGenerator.hpp"
#include "grid/GridStorage.hpp"

#include "sgpp.hpp"

namespace sg
{

StretchedTrapezoidBoundaryGridGenerator::StretchedTrapezoidBoundaryGridGenerator(GridStorage* storage) : storage(storage)
{
}

StretchedTrapezoidBoundaryGridGenerator::~StretchedTrapezoidBoundaryGridGenerator()
{
}

void StretchedTrapezoidBoundaryGridGenerator::regular(size_t level)
{
	HashGenerator gen;
	gen.regularWithBoundaries(this->storage, level, true);
}

void StretchedTrapezoidBoundaryGridGenerator::refine(RefinementFunctor* func)
{
	HashRefinementBoundaries refine;
	refine.free_refine(this->storage, func);
}

size_t StretchedTrapezoidBoundaryGridGenerator::getNumberOfRefinablePoints()
{
	HashRefinementBoundaries refine;
	return refine.getNumberOfRefinablePoints(this->storage);
}

void StretchedTrapezoidBoundaryGridGenerator::coarsen(CoarseningFunctor* func, DataVector* alpha)
{
	HashCoarsening coarsen;
	coarsen.free_coarsen(this->storage, func, alpha);
}

void StretchedTrapezoidBoundaryGridGenerator::coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly)
{
	HashCoarsening coarsen;
	coarsen.free_coarsen_NFirstOnly(this->storage, func, alpha, numFirstOnly);
}

size_t StretchedTrapezoidBoundaryGridGenerator::getNumberOfRemoveablePoints()
{
	HashCoarsening coarsen;
	return coarsen.getNumberOfRemovablePoints(this->storage);
}

void StretchedTrapezoidBoundaryGridGenerator::refineMaxLevel(RefinementFunctor* func, unsigned int maxLevel)
{
	HashRefinementBoundariesMaxLevel refine;
	refine.refineToMaxLevel(this->storage, func, maxLevel);
}

size_t StretchedTrapezoidBoundaryGridGenerator::getNumberOfRefinablePointsToMaxLevel(unsigned int maxLevel)
{
	HashRefinementBoundariesMaxLevel refine;
	return refine.getNumberOfRefinablePointsToMaxLevel(this->storage, maxLevel);
}

}
