/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

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
