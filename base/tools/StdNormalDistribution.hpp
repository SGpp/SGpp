/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef STDNORMALDISTRIBUTION_HPP
#define STDNORMALDISTRIBUTION_HPP

#include <cmath>

namespace sg {
  namespace base {

    /**
     * This provides an implementation of the standradized normal
     * distribution.
     *
     * Also a density function for the normal distribution is provided!
     */
    class StdNormalDistribution {
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

        /*
         * Calculates the Density values of the normal distribution
         *
         * taken from http://www.richelbilderbeek.nl/CppGetDensityNormal.htm
         *
         * @param x the value for which the density value should be calculated
         * @param mu the expected value of the normal distribution
         * @param sigma the standard deviation of the normal distribution
         */
        double getDensity(const double x, const double mu, const double sigma);


        /*
         * Calculates the Density values of the normal distribution (normed variant)
         *
         * @param x the value for which the density value should be calculated
         * @param mu the expected value of the normal distribution
         * @param sigma the standard deviation of the normal distribution
         */
        double getNormedDensity(const double x, const double mu, const double sigma);
    };

  }
}

#endif /* STDNORMALDISTRIBUTION_HPP */
