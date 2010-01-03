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

#include "tools/common/StdNormalDistribution.hpp"

namespace sg
{

StdNormalDistribution::StdNormalDistribution()
{
}

StdNormalDistribution::~StdNormalDistribution()
{
}

double StdNormalDistribution::getCumulativeDensity(const double x)
{
	const double c0 = 0.2316419;
	const double c1 = 1.330274429;
	const double c2 = 1.821255978;
	const double c3 = 1.781477937;
	const double c4 = 0.356563782;
	const double c5 = 0.319381530;
	const double c6 = 0.398942280401;

	const double negative = (x < 0 ? 1.0 : 0.0);
	const double xPos = (x < 0.0 ? -x : x);
	const double k = 1.0 / ( 1.0 + (c0 * xPos));
	const double y1 =(((((((c1*k-c2)*k)+c3)*k)-c4)*k)+c5)*k;
	const double y2 = 1.0 - (c6*std::exp(-0.5*xPos*xPos)*y1);

	return ((1.0-negative)*y2) + (negative*(1.0-y2));
}


double StdNormalDistribution::getDensity(const double x)
{
	const double mean = 0.0;
	const double stddev = 1.0;

	const double firstTerm = 1.0 / (stddev * std::sqrt(2.0 * M_PI));
	const double secondTerm = -( (x - mean) * (x - mean) / (2.0 * stddev * stddev) );
	const double result = firstTerm * std::exp(secondTerm);

	return result;
}

}
