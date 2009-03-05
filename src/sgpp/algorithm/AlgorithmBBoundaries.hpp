/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
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
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef ALGORTIHMBBOUNDARIES_HPP
#define ALGORTIHMBBOUNDARIES_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

#include "algorithm/GetAffectedBasisFunctionsBoundaries.hpp"

#include <vector>
#include <utility>
#include <iostream>

#ifdef USEOMP
#include <omp.h>
#endif

namespace sg {

/**
 * Basic multiplaction with B and B^T.
 *
 * Boundaries are supported
 */
template<class BASIS>
class AlgorithmBBoundaries
{
public:

	void mult(GridStorage* storage, BASIS& basis, DataVector& source, DataVector& x, DataVector& result)
	{
		typedef std::vector<std::pair<size_t, double> > IndexValVector;

		result.setAll(0.0);
		size_t source_size = source.getSize();

		std::vector<double> line;
		IndexValVector vec;

		GetAffectedBasisFunctionsBoundaries<BASIS> ga(storage);

		for(size_t i = 0; i < source_size; i++)
		{
			vec.clear();

			x.getLine(i, line);

			ga(basis, line, vec);

			for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
			{
				result[iter->first] += iter->second * source[i];
			}
		}
	}

	void mult_transpose(GridStorage* storage, BASIS& basis, DataVector& source, DataVector& x, DataVector& result)
	{
		typedef std::vector<std::pair<size_t, double> > IndexValVector;

		result.setAll(0.0);
		size_t result_size = result.getSize();

		std::vector<double> line;
		IndexValVector vec;

		GetAffectedBasisFunctionsBoundaries<BASIS> ga(storage);

#ifdef USEOMP
		#pragma omp parallel for private (vec, line) shared (result) schedule (static)
		for(size_t i = 0; i < result_size; i++)
		{
			vec.clear();

			x.getLine(i, line);

			ga(basis, line, vec);

			for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
			{
				result[i] += iter->second * source[iter->first];
			}
		}
#else
		for(size_t i = 0; i < result_size; i++)
		{
			vec.clear();

			x.getLine(i, line);

			ga(basis, line, vec);

			for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
			{
				result[i] += iter->second * source[iter->first];
			}
		}
#endif /* USEOMP */
	}


protected:

};

}

#endif /* ALGORTIHMBBOUNDARIES_HPP */
