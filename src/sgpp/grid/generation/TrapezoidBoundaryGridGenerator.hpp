/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
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

#ifndef TRAPEZOIDBOUNDARYGRIDGENERATOR_HPP
#define TRAPEZOIDBOUNDARYGRIDGENERATOR_HPP

#include "grid/GridStorage.hpp"
#include "grid/generation/GridGenerator.hpp"

namespace sg
{

/**
 * This class provides the interface for the grid generation
 * for grids with boundaries, pentagon cut through sub space scheme
 */
class TrapezoidBoundaryGridGenerator : public GridGenerator
{
public:
	/**
	 * Constructor
	 *
	 * @param storage template type that holds the grid points
	 */
	TrapezoidBoundaryGridGenerator(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~TrapezoidBoundaryGridGenerator();

	virtual void regular(size_t level);
	virtual void refine(RefinementFunctor* func);
	virtual int getNumberOfRefinablePoints();

protected:
	/// Pointer to the grid's storage object
	GridStorage* storage;
};

}

#endif /* TRAPEZOIDBOUNDARYGRIDGENERATOR_HPP */
