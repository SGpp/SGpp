/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/******************************************************************************/

#include "grid/generation/BoundaryUScaledGridGenerator.hpp"
#include "grid/GridStorage.hpp"

#include "sgpp.hpp"

namespace sg
{

/**
 * Constructor
 *
 * @param storage template type that holds the grid points
 */
BoundaryUScaledGridGenerator::BoundaryUScaledGridGenerator(GridStorage* storage) : storage(storage)
{
}

/**
 * Destructor
 */
BoundaryUScaledGridGenerator::~BoundaryUScaledGridGenerator()
{
}

/**
 * creates a regular grid with boundaries, pentagon cut
 *
 * @param level maximum level of the grid
 */
void BoundaryUScaledGridGenerator::regular(size_t level)
{
	HashGenerator gen;
	gen.regularWithBoundaries(this->storage, level, true);
}

/**
 * refines a regular grid
 *
 * @param func pointer to refinement function
 */
void BoundaryUScaledGridGenerator::refine(RefinementFunctor* func)
{
	HashRefinement refine;
	refine.free_refine(this->storage, func, true);
}

}
