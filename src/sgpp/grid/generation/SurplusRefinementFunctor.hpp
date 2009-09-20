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

#ifndef SURPLUSREFINEMENTFUNCTOR_HPP
#define SURPLUSREFINEMENTFUNCTOR_HPP

#include "data/DataVector.hpp"
#include "grid/generation/RefinementFunctor.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * This abstracts the refinement criteria out of the refinement algorithm
 * @version $HEAD$ 
 */
class SurplusRefinementFunctor : public RefinementFunctor
{
public:
	/**
	 * Constructor
	 *
	 * @param alpha DataVector that is basis for refinement decisions
	 * @param refinements_num number of grid points which should be refined
	 * @param threshold
	 */
	SurplusRefinementFunctor(DataVector* alpha, int refinements_num = 1, double threshold = 0.0) : alpha(alpha), refinements_num(refinements_num), threshold(threshold)
	{
	}

	/**
	 * Destructor
	 */
	virtual ~SurplusRefinementFunctor() {}


	virtual double operator()(GridStorage* storage, size_t seq)
	{
		return fabs(alpha->get(seq));
	}

	virtual double start()
	{
		return 0.0;
	}

	int getRefinementsNum()
	{
		return this->refinements_num;
	}
	
	double getRefinementThreshold()
	{
		return this->threshold;
	}


protected:
	/// pointer to the vector that stores the alpha values
	DataVector* alpha;
	
	/// number of grid points to refine
	int refinements_num;
	
	/// threshold, only the points with greater to equal absolute values of the refinement criterion (e.g. alpha or error) will be refined
	double threshold;
};

}

#endif /* SURPLUSREFINEMENTFUNCTOR_HPP */
