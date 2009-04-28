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

#ifndef REFINEMENTFUNCTOR_HPP
#define REFINEMENTFUNCTOR_HPP

#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * Abstract class the defines the interface that refinement functors have to provide
 */
class RefinementFunctor
{
public:
	typedef double value_type;

	/**
	 * Constructor
	 */
	RefinementFunctor() {}

	/**
	 * Destructor
	 */
	virtual ~RefinementFunctor() {}

	/**
	 * This should be returning a refinement value for every grid point.
	 * The point with the highest value will be refined.
	 *
	 * @param storage pointer to the grids storage object
	 * @param seq sequence number in the coefficients array
	 *
	 * @return refinement value
	 */
	virtual double operator()(GridStorage* storage, size_t seq) = 0;

	/**
	 * This should return a bottom.
	 *
	 * @return the minimal value
	 */
	virtual double start() = 0;

	/**
	 * makes it possible to refine on several points
	 *
	 * @return number of refinements ????
	 */
	virtual int getRefinementsNum(){ return 1;}
};

}

#endif /* REFINEMENTFUNCTOR_HPP */
