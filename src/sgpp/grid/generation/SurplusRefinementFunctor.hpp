/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SURPLUSREFINEMENTFUNCTOR_HPP
#define SURPLUSREFINEMENTFUNCTOR_HPP

#include "data/DataVector.hpp"
#include "grid/generation/RefinementFunctor.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * A refinement functor, refining according to the maximal absolute values in a DataVector provided.
 * @version $HEAD$ 
 */
class SurplusRefinementFunctor : public RefinementFunctor
{
public:
	/**
	 * Constructor.
	 *
	 * @param alpha DataVector that is basis for refinement decisions. The i-th entry corresponds to the i-th grid point.
	 * @param refinements_num Number of grid points which should be refined (if possible - there could be less refinable grid points)
	 * @param threshold The absolute value of the entries have to be greater or equal than the threshold
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
