/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (dirk.pflueger@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef GRIDGENERATOR_HPP
#define GRIDGENERATOR_HPP

#include "grid/generation/RefinementFunctor.hpp"

namespace sg
{

/**
 * Abstract class that defines the interfaces for the different grid's GridGenerators
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
	 * Creates a regular grid for a certain level @f$ n @f$, i.e., @f$ V_n^{(1)} = \bigoplus_{|\vec{l}|_1 \leq n+d-1} W_{\vec{l}}@f$.
	 *
	 * @param level maximum level of the grid
	 */
	virtual void regular(size_t level) = 0;

	/**
	 * Refines a regular grid according to the settings of the RefinementFunctor func.
	 *
	 * @param func pointer to refinement functor
	 */
	virtual void refine(RefinementFunctor* func) = 0;

	/**
	 * Returns the number of points on the grid that can be refined in the next iteration
	 * 
	 * @return the number of points on the grid that can be refined
	 */
	virtual int getNumberOfRefinablePoints() = 0;
};

}

#endif /* GRIDGENERATOR_HPP */
