/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/generation/StandardGridGenerator.hpp"
#include "grid/GridStorage.hpp"

#include "exception/generation_exception.hpp"

#include "sgpp.hpp"

namespace sg
{

StandardGridGenerator::StandardGridGenerator(GridStorage* storage) : storage(storage)
{
}

StandardGridGenerator::~StandardGridGenerator()
{
}

void StandardGridGenerator::regular(size_t level)
{
	HashGenerator gen;
	gen.regular(this->storage, level);
}

void StandardGridGenerator::refine(RefinementFunctor* func)
{
	HashRefinement refine;
	refine.free_refine(this->storage, func);
}

int StandardGridGenerator::getNumberOfRefinablePoints()
{
	HashRefinement refine;
	return refine.getNumberOfRefinablePoints(this->storage);
}

}
