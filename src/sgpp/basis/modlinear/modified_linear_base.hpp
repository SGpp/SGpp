/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (dirk.pflueger@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef MODIFIED_LINEAR_BASE_HPP
#define MODIFIED_LINEAR_BASE_HPP

#include <cmath>

namespace sg
{

/**
 * modified linear base functions.
 *
 * @version $HEAD$
 */
 template<class LT, class IT>
class modified_linear_base
{
public:
	/**
	 * Evaluate a basis function.
	 * Has a dependence on the absolute position of grid point and support.
	 */
	double eval(LT level, IT index, double p)
	{
		if(level == 1)
		{
			return 1.0;
		}
		else if(index == 1)
		{
			return 2.0 - (1<<level) * p;
		}
		else if(static_cast<int>(index) == static_cast<int>((1<<level)-1))
		{
			return (1<<level) * p - index + 1.0;
		}
		return 1.0 - fabs((1<<level) * p - index);
	}
};

}

#endif /* MODIFIED_LINEAR_BASE_HPP */

