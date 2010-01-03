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

#ifndef STDNORMALDISTRIBUTION_HPP
#define STDNORMALDISTRIBUTION_HPP

#include <cmath>

namespace sg
{

class StdNormalDistribution
{
public:
	/**
	 * Std-Constructor
	 */
	StdNormalDistribution();

	/**
	 * Std-Destructor
	 */
	~StdNormalDistribution();

	/**
	 * Calculates the Cumulative Density values of the standard normal distribution
	 * (expected values = 0.0, standard deviation = 1.0)
	 *
	 * taken from http://www.richelbilderbeek.nl/CppGetCumulativeDensityNormal.htm
	 *
	 * @param x the value for which the cumulative density value should be calculated
	 */
	double getCumulativeDensity(const double x);

	/*
	 * Calculates the Density values of the standard normal distribution
	 * (expected values = 0.0, standard deviation = 1.0)
	 *
	 * taken from http://www.richelbilderbeek.nl/CppGetDensityNormal.htm
	 *
	 * @param x the value for which the density value should be calculated
	 */
	double getDensity(const double x);
};

}

#endif /* STDNORMALDISTRIBUTION_HPP */
