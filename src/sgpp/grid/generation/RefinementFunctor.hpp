/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008-2009 Dirk Pflueger (pflueged@in.tum.de)                */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
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
 * @version $HEAD$  
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
	 * The point with the highest value will be refined first.
	 *
	 * @param storage pointer to the grids storage object
	 * @param seq sequence number in the coefficients array
	 *
	 * @return refinement value
	 */
	virtual double operator()(GridStorage* storage, size_t seq) = 0;

	/**
	 * This should return the initial value of refinement criterion (e.g. alpha or error).
	 *
	 * @return the initial value
	 */
	virtual double start() = 0;

	/**
	 * Returns the maximal number of points that should be refined.
	 * 
	 * The maximal number of points to refine is set in the constructor of implementation class. 
	 *
	 * @return number of points that should refined. Default value: 1.
	 */
	virtual int getRefinementsNum(){ return 1;}
	
	/**
	 * Returns the threshold value.
	 * 
	 * Only the grid points with absolute value of refinement criterion (e.g. alpha or error) greater 
	 * or equal to this threshold will be refined
	 * 
	 * @return threshold value for refinement. Default value: 0.
	 */
	virtual double getRefinementThreshold() = 0;
};

}

#endif /* REFINEMENTFUNCTOR_HPP */
