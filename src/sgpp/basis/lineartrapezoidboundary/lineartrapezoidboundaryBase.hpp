/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
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

#ifndef LINEARTRAPEZOIDBOUNDARYBASE_HPP
#define LINEARTRAPEZOIDBOUNDARYBASE_HPP

#include <cmath>

namespace sg
{

/**
 * linear basis functions with boundaries
 * And here we have another implicit dependence on tensor products
 *
 * @version $HEAD$
 */
template<class LT, class IT>
class lineartrapezoidboundaryBase
{
public:
	/**
	 * Evaluate a basis function.
	 * Has a dependence on the absolute position of grid point and support.
	 *
	 * @param level the level of the current basis function
	 * @param index the index of the current basis function
	 * @param p the absolute position of the evaluation point
	 */
	double eval(LT level, IT index, double p)
	{
		if (level == 0)
		{
			if (index == 0)
			{
				return 1.0 - p;
			}
			if (index == 1)
			{
				return p;
			}
		}
		else
		{
			return 1.0 - fabs((1<<level) * p - index);
		}
		// should not happen
		return 0.0;
	}

	/**
	 * Evaluate a basis function with an offset and scaling factor
	 * Has a dependence on the absolute position of grid point and support.
	 *
	 * @param level the level of the current basis function
	 * @param index the index of the current basis function
	 * @param p the absolute position of the evaluation point
	 * @param q the scaling factor of the basis function
	 * @param t the offset of the basis function
	 */
	double eval(LT level, IT index, double p, double q, double t)
	{
		if (level == 0)
		{
			if (index == 0)
			{
				return 1.0 - ((1.0/q)*(p-t));
			}
			if (index == 1)
			{
				return ((1.0/q)*(p-t));
			}
		}
		else
		{
			return 1.0 - ((1.0/q)*(fabs(((1<<level)*(p-t))-(q*index))));
		}
		// should not happen
		return 0.0;
	}
};

}

#endif /* LINEARTRAPEZOIDBOUNDARYBASE_HPP */
