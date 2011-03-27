/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/generation/BoundaryGridGenerator.hpp"
#include "grid/GridStorage.hpp"

#include "sgpp.hpp"

namespace sg
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

int BoundaryGridGenerator::getNumberOfRefinablePoints()
{
	HashRefinementBoundaries refine;
	return refine.getNumberOfRefinablePoints(this->storage);
}

}
