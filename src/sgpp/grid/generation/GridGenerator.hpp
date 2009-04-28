/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef GRIDGENERATOR_HPP
#define GRIDGENERATOR_HPP

#include "grid/generation/RefinementFunctor.hpp"

namespace sg
{

/**
 * abstract class that defines the interfaces for the different grid's GridGenerators
 */
class GridGenerator
{
public:
	/**
	 * Constructor
	 */
	GridGenerator() {}

	/**
	 * Destructor
	 */
	virtual ~GridGenerator() {}

	/**
	 * creates a regular grid
	 *
	 * @param level maximum level of the grid
	 */
	virtual void regular(size_t level) = 0;

	/**
	 * refines a regular grid
	 *
	 * @param func pointer to refinement function
	 */
	virtual void refine(RefinementFunctor* func) = 0;

	/**
	 * @todo add comment
	 */
	virtual int getNumberOfRefinablePoints() = 0;
};

}

#endif /* GRIDGENERATOR_HPP */
