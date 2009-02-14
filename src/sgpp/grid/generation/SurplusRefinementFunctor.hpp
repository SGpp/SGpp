/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
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
/*****************************************************************************/

#ifndef SURPLUSREFINEMENTFUNCTOR_HPP
#define SURPLUSREFINEMENTFUNCTOR_HPP

#include "data/DataVector.h"
#include "grid/generation/RefinementFunctor.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * This abstracts the refinement criteria out of the refinement algorithm
 */
class SurplusRefinementFunctor : public RefinementFunctor
{
public:


	SurplusRefinementFunctor(DataVector* alpha) : alpha(alpha)
	{
	}

	virtual ~SurplusRefinementFunctor() {}

	/**
	 * This should be returning a refinement value for every grid point.
	 * The point with the highest value will be refined.
	 */
	virtual double operator()(GridStorage* storage, size_t seq)
	{
		return fabs(alpha->get(seq));
	}

	/**
	 * This should return a bottom.
	 */
	virtual double start()
	{
		return 0.0;
	}


protected:
	DataVector* alpha;
};

}

#endif /* SURPLUSREFINEMENTFUNCTOR_HPP */
