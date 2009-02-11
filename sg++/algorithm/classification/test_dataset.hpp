/*****************************************************************************/
/* This file is part of sg++, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sg++ is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sg++ is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sg++; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef TEST_DATASET_HPP
#define TEST_DATASET_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

#include <vector>
#include <utility>
#include <iostream>

namespace sg {

/**
 * Returns the number of correctly classified instances in data
 */
template<class BASIS>
double test_dataset( GridStorage* storage, BASIS& basis, DataVector& alpha, DataVector& data, DataVector& classes)
{
	typedef std::vector<std::pair<size_t, double> > IndexValVector;

	double correct = 0;

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

	return correct;

}

}

#endif /* TEST_DATASET_HPP */
