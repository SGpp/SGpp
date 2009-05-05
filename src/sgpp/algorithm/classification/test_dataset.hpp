/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
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

#ifndef TEST_DATASET_HPP
#define TEST_DATASET_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#include <vector>
#include <utility>
#include <iostream>

#ifdef USEOMP
#include <omp.h>
#endif

namespace sg {

/**
 * Returns the number of correctly classified instances in data without boundaries
 *
 * @param storage GridStorage object that contains the grid points
 * @param basis reference to class that implements to current basis
 * @param alpha the coefficients of the grid points
 * @param data the data the should be tested
 * @param classes the classes computed by the sparse grid's classification algorithm
 */
template<class BASIS>
double test_dataset( GridStorage* storage, BASIS& basis, DataVector& alpha, DataVector& data, DataVector& classes)
{
	typedef std::vector<std::pair<size_t, double> > IndexValVector;

	double correct = 0;

#ifdef USEOMP
	#pragma omp parallel shared(correct)
	{
		size_t size = data.getSize();

		std::vector<double> point;

		GetAffectedBasisFunctions<BASIS> ga(storage);

		#pragma omp for schedule(static)
		for(size_t i = 0; i < size; i++)
		{

			IndexValVector vec;
			double result = 0;

			data.getLine(i, point);

			ga(basis, point, vec);

			for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
			{
				result += iter->second * alpha[iter->first];
			}

			#pragma omp critical
			{
				if( (result >= 0 && classes[i] >= 0) || (result < 0 && classes[i] < 0) )
				{
					correct++;
				}
			}
		}
	}
#else
	size_t size = data.getSize();

	std::vector<double> point;

	GetAffectedBasisFunctions<BASIS> ga(storage);

	for(size_t i = 0; i < size; i++)
	{

		IndexValVector vec;
		double result = 0;

		data.getLine(i, point);

		ga(basis, point, vec);

		for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
		{
			result += iter->second * alpha[iter->first];
		}

		if( (result >= 0 && classes[i] >= 0) || (result < 0 && classes[i] < 0) )
		{
			correct++;
		}
	}
#endif

	return correct;

}

}

#endif /* TEST_DATASET_HPP */
