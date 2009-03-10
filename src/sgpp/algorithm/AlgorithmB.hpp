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

#ifndef ALGORTIHMB_HPP
#define ALGORTIHMB_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

#include "algorithm/GetAffectedBasisFunctions.hpp"

#include <vector>
#include <utility>
#include <iostream>

#ifdef USEOMP
#include <omp.h>
#endif

namespace sg {

/**
 * Basic multiplaction with B and B^T.
 * This is just wrong ... there should be some functor for the BASIS type
 */
template<class BASIS>
class AlgorithmB
{
public:

	void mult(GridStorage* storage, BASIS& basis, DataVector& source, DataVector& x, DataVector& result)
	{
		typedef std::vector<std::pair<size_t, double> > IndexValVector;

		result.setAll(0.0);
#ifdef USEOMP
		#pragma omp parallel shared(result)
		{
			size_t source_size = source.getSize();

			std::vector<double> line;
			IndexValVector vec;

			GetAffectedBasisFunctions<BASIS> ga(storage);

			#pragma omp for schedule(static)
			for(size_t i = 0; i < source_size; i++)
			{
				vec.clear();

				x.getLine(i, line);

				ga(basis, line, vec);

				#pragma omp critical
				{
					for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
					{
						result[iter->first] += iter->second * source[i];
					}
				}
			}
		}
#else
		size_t source_size = source.getSize();

		std::vector<double> line;
		IndexValVector vec;

		GetAffectedBasisFunctions<BASIS> ga(storage);

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
#endif /* USEOMP */
	}

	void mult_transpose(GridStorage* storage, BASIS& basis, DataVector& source, DataVector& x, DataVector& result)
	{
		typedef std::vector<std::pair<size_t, double> > IndexValVector;

		result.setAll(0.0);

#ifdef USEOMP
		#pragma omp parallel shared (result)
		{
			size_t result_size = result.getSize();

			std::vector<double> line;
			IndexValVector vec;

			GetAffectedBasisFunctions<BASIS> ga(storage);

			#pragma omp for schedule (static)
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
		}
#else
		size_t result_size = result.getSize();

		std::vector<double> line;
		IndexValVector vec;

		GetAffectedBasisFunctions<BASIS> ga(storage);

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

#endif /* ALGORTIHMB_HPP */
