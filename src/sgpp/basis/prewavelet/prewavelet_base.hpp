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
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef PREWAVELET_BASE_HPP
#define PREWAVELET_BASE_HPP

#include <cmath>
#include <iostream>

namespace sg
{
namespace base
{

/**
 * Class representing a prewavelet base.
 * A prewavelet \f$\psi\f$ is a combination of 5 normal hat functions \f$\phi\f$
 * with a \f$\left[\frac{1}{10}, -\frac{6}{10}, 1, -\frac{6}{10}, \frac{1}{10}\right]\f$ stamp:
 * \f[
 * \psi_{(l,i)} = \frac{1}{10}\phi_{(l,i-2)} - \frac{6}{10}\phi_{(l,i-1)} + \phi_{(l,i)} - \frac{6}{10} \phi_{(l,i+1)} + \frac{1}{10}\phi_{(l,i+2)}
 * \f]
 * Next to the border:
 * \f[
 * \psi_{(l,1)} = \frac{9}{10}\phi_{(l,1)} - \frac{6}{10}\phi_{(l,2)} + \frac{1}{10}\phi_{(l,3)}
 * \f]
 * For level \f$ l = 1\f$ the prewavelet and linear basis are equivalent.
 * \image html prewavelets.png "Normal prewavelet base function and a left border prewavelet. The bold dots indicate the position of other prewavelet basis functions on the same level, the thin dots representing grid points missing in the sparse grid on that level. Please note that the left and right TWO neighbors are interfering with a specific prewavelet base."
 * The prewavelets form a semi-orthogonal basis. That means
 * \f$<\psi_{(i,l)},\psi_{(j,k)}> = 0\f$ if \f$l\neq k\f$. On the same level, the prewavelets are not
 * orthogonal. This property will ease the calculation of some specific operations, but on the other hand,
 * this advantage is bought with a wider support of the ansatzfunctions.
 */
template<class LT, class IT>
class prewavelet_base
{
private:

	static const double normal_stamp[];
	static const double border_stamp[];

	double inline evalNormalHat(LT level, IT index, double p)
	{
		return 1.0 - fabs((1 << level) * p - index);
	}

public:

	/**
	 * Evaluate a base functions.
	 * Has a dependence on the absolute position of grid point and support.
	 */
	double eval(LT level, IT index, double p)
	{
		/* A prewavelet is simply the sum of neighboring hat functions with a certain
		 * stamp. So we have to figure out, which of the 5 hat functions is the one
		 * "hit" by p, modify the index and use the normal hat-function calculation.
		 */

//		//First, check, if the given point is really in the support area of this basis-function!
//		int int_index = index; // needed to avoid overrun with index-3
//		if ((1.0 / (1 << level)) * (int_index - 3) > p || (1.0 / (1 << level))
//				* (index + 3) < p)
//		{
//			return 0.0;
//		}

		if (p == 0.0 || p == 1.0)
		{
			return 0.0;
		}

		if (1 == level)
		{
			return evalNormalHat(level, index, p);
		}
		else if (1 == index) //left border
		{
			//Index of the affected hatbasis. The affected bases are ab and ab + 1
			int ab = floor(p * (1 << level));
			if (ab == 0)
			{
				return 0.9 * evalNormalHat(level, 1, p);
			}
			else if (ab == 3)
			{
				return 0.1 * evalNormalHat(level, 3, p);
			}
			else
			{
				return border_stamp[ab - 1] * evalNormalHat(level, ab, p)
						+ border_stamp[ab] * evalNormalHat(level, ab + 1, p);
			}
		}
		else if ((unsigned int)(1 << level) - 1 == index) //right border
		{
			//Index of the affected hatbasis. The affected bases are ab and ab + 1
			int ab = floor(p * (1 << level));
			if (ab == (1 << level) - 1)
			{
				return 0.9 * evalNormalHat(level, (1 << level) - 1, p);
			}
			else if (ab == (1 << level) - 4)
			{
				return 0.1 * evalNormalHat(level, (1 << level) - 3, p);
			}
			else
			{
				return border_stamp[(1 << level) - 1 - ab] * evalNormalHat(
						level, ab, p) + border_stamp[(1 << level) - 2 - ab]
						* evalNormalHat(level, ab + 1, p);
			}
		}
		else
		{
			//Index of the affected hatbasis. The affected bases are ab and ab + 1
			unsigned int ab = floor(p * (1 << level));
			if (ab == index - 3)
			{
				return normal_stamp[0] * evalNormalHat(level, ab + 1, p);
			}
			else if (ab == index + 2)
			{
				return normal_stamp[4] * evalNormalHat(level, ab, p);
			}
			else
			{
				int stamp = ab - index + 2;
				return normal_stamp[stamp] * evalNormalHat(level, ab, p)
						+ normal_stamp[stamp + 1] * evalNormalHat(level,
								ab + 1, p);
			}
		}
	}
};

template<class LT, class IT>
const double prewavelet_base<LT, IT>::normal_stamp[] =
{ 0.1, -0.6, 1.0, -0.6, 0.1 };

template<class LT, class IT>
const double prewavelet_base<LT, IT>::border_stamp[] =
{ 0.9, -0.6, 0.1 };

}
}

#endif /* PREWAVELET_BASE_HPP */
