/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SURPLUSCOARSENINGFUNCTOR_HPP
#define SURPLUSCOARSENINGFUNCTOR_HPP

#include "data/DataVector.hpp"
#include "grid/generation/CoarseningFunctor.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * A coarsening functor, removing points according to the minimal absolute values in a DataVector provided.
 * @version $HEAD$
 */
class SurplusCoarseningFunctor : public CoarseningFunctor
{
public:
	/**
	 * Constructor.
	 *
	 * @param alpha DataVector that is basis for coarsening decisions. The i-th entry corresponds to the i-th grid point.
	 * @param removements_num Number of grid points which should be removed (if possible - there could be less removeabled grid points)
	 * @param threshold The absolute value of the entries have to be greater or equal than the threshold
	 */
	SurplusCoarseningFunctor(DataVector* alpha, int removements_num = 1, double threshold = 0.0) : alpha(alpha), removements_num(removements_num), threshold(threshold)
	{
	}

	/**
	 * Destructor
	 */
	virtual ~SurplusCoarseningFunctor() {}


	virtual double operator()(GridStorage* storage, size_t seq)
	{
		return fabs(alpha->get(seq));
	}

	virtual double start()
	{
		return 0.0;
	}

	int getRemovementsNum()
	{
		return this->removements_num;
	}

	double getCoarseningThreshold()
	{
		return this->threshold;
	}


protected:
	/// pointer to the vector that stores the alpha values
	DataVector* alpha;

	/// number of grid points to remove
	int removements_num;

	/// threshold, only the points with greater to equal absolute values of the refinement criterion (e.g. alpha or error) will be refined
	double threshold;
};

}

#endif /* SURPLUSCOARSENINGFUNCTOR_HPP */
