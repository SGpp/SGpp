/* ****************************************************************************
 * Copyright (C) 2014 Universitaet Stuttgart                                   *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Andreas Doerr, Marcel Schneider, Matthias Moegerle
#ifndef SOBOLSAMPLEGENERATOR_H_
#define SOBOLSAMPLEGENERATOR_H_

#include "base/datatypes/DataVector.hpp"
#include "mcm/SampleGenerator.hpp"

namespace sg {
namespace mcm {

/**
 * The class SobolSsampleGenerator implements the Sobol Sequence.
 * Primitive Polynomials are taken from http://web.maths.unsw.edu.au/~fkuo/sobol/
 * The output is based on a unit cube, each dimension from 0.0 to 1.0.
 */
class SobolSampleGenerator: public SampleGenerator {

public:

	SobolSampleGenerator(size_t dimensions);
	/**
	 * @brief Constructor, init
	 *
	 * @param dimensions define Dimensions
	 * @param seed seed skips the first x values in sobol' sequence
	 *
	 * @return void
	 */
	SobolSampleGenerator(size_t dimensions, size_t seed);

	virtual ~SobolSampleGenerator() {};

	/**
	 *
	 *  Example:
	 *
	 *       N    Binary     BIT
	 *    ----    --------  ----
	 *       0           0     0
	 *       1           1     1
	 *       2          10     2
	 *       3          11     2
	 *       4         100     3
	 *       5         101     3
	 *       6         110     3
	 *       7         111     3
	 *       8        1000     4
	 *       9        1001     4
	 *      10        1010     4
	 *      11        1011     4
	 *      12        1100     4
	 *      13        1101     4
	 *      14        1110     4
	 *      15        1111     4
	 *      16       10000     5
	 *      17       10001     5
	 *    1023  1111111111    10
	 *    1024 10000000000    11
	 *    1025 10000000001    11
	 *
	 *      Licensing:
	 *
	 *            This code is distributed under the GNU LGPL license.
	 *
	 *      Modified:
	 *
	 *            22 February 2011
	 *
	 *    Author:
	 *
	 *        Original MATLAB version by John Burkardt.
	 *        PYTHON version by Corrado Chisari
	 *        C++ version by Andreas Doerr, Marcel Schneider, Matthias Moegerle
	 *
	 *      Parameters:
	 *	     @brief I4_BIT_HI1 returns the position of the high 1 bit base 2 in an integer.
	 *       @param n the integer to be measured.
	 *               N should be nonnegative.  If N is nonpositive, the value will always be 0.
	 *	     @return Output, integer BIT, the number of bits base 2.
	 *
	 *
	 *
	 */
	size_t i4_bit_hi1(size_t n);

	/*****************************************************************************80
	 *
	 *  Example:
	 *
	 *       N    Binary     BIT
	 *    ----    --------  ----
	 *       0           0     1
	 *       1           1     2
	 *       2          10     1
	 *       3          11     3
	 *       4         100     1
	 *       5         101     2
	 *       6         110     1
	 *       7         111     4
	 *       8        1000     1
	 *       9        1001     2
	 *      10        1010     1
	 *      11        1011     3
	 *      12        1100     1
	 *      13        1101     2
	 *      14        1110     1
	 *      15        1111     5
	 *      16       10000     1
	 *      17       10001     2
	 *    1023  1111111111     1
	 *    1024 10000000000     1
	 *    1025 10000000001     1
	 *
	 *      Licensing:
	 *
	 *    This code is distributed under the GNU LGPL license.
	 *
	 *      Modified:
	 *
	 *            22 February 2011
	 *
	 *    Author:
	 *
	 *        Original MATLAB version by John Burkardt.
	 *        PYTHON version by Corrado Chisari
	 *        C++ version by Andreas Doerr, Marcel Schneider, Matthias Moegerle
	 *
	 *  Parameters:
	 *	     @brief I4_BIT_LO0 returns the position of the low 0 bit base 2 in an integer.
	 *        @param Input, integer N, the integer to be measured.
	 *            N should be nonnegative.
	 *	     @return size_t BIT, the position of the low 1 bit.
	 */
	size_t i4_bit_lo0(size_t n);

	/*****************************************************************************80
	 *
	 *    Discussion:
	 *
	 *        The routine adapts the ideas of Antonov and Saleev->
	 *
	 *    Licensing:
	 *
	 *        This code is distributed under the GNU LGPL license.
	 *
	 *    Modified:
	 *
	 *            22 February 2011
	 *
	 *    Author:
	 *
	 *        Original FORTRAN77 version by Bennett Fox.
	 *        MATLAB version by John Burkardt.
	 *        PYTHON version by Corrado Chisari
	 *        C++ version by Andreas Doerr, Marcel Schneider, Matthias Moegerle
	 *
	 *    Reference:
	 *
	 *        Antonov, Saleev,
	 *        USSR Computational Mathematics and Mathematical Physics,
	 *        Volume 19, 1980, pages 252 - 256.
	 *
	 *        Paul Bratley, Bennett Fox,
	 *        Algorithm 659:
	 *        Implementing Sobol's Quasirandom Sequence Generator,
	 *        ACM Transactions on Mathematical Software,
	 *        Volume 14, Number 1, pages 88-100, 1988.
	 *
	 *        Bennett Fox,
	 *        Algorithm 647:
	 *        Implementation and Relative Efficiency of Quasirandom
	 *        Sequence Generators,
	 *        ACM Transactions on Mathematical Software,
	 *        Volume 12, Number 4, pages 362-376, 1986.
	 *
	 *        Ilya Sobol,
	 *        USSR Computational Mathematics and Mathematical Physics,
	 *        Volume 16, pages 236-242, 1977.
	 *
	 *        Ilya Sobol, Levitan,
	 *        The Production of Points Uniformly Distributed in a Multidimensional
	 *        Cube (in Russian),
	 *        Preprint IPM Akad. Nauk SSSR,
	 *        Number 40, Moscow 1976.
	 *
	 *    Parameters:
	 *
	 *	     @brief getSample generates a new quasirandom Sobol vector with each call.
	 *        @param dv Vector which in which results will be written
	 *	     @return void.
	 */
	void getSample(sg::base::DataVector& sample);

private:

	// seed skip the first x samples in sobols' sequence
	size_t seed;

	// index number of current sample [1..numberOfSamples]
	size_t numberOfCurrentSample;

	/**
	 * @brief Method returns the total number of samples which can be generated
	 * according to the sample generator settings (dimensions and subdivision into strata)
	 *
	 * @return size_t total number of samples
	 */
	double** sobol_points(char *dir_file);

};

}
}

#endif /* SOBOLSAMPLEGENERATOR_H_ */
