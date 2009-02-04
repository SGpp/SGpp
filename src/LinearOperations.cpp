/*
This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems

Copyright (C) 2008 Joerg Blank (blankj@in.tum.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "linear/LinearOperations.hpp"

#include "sgpp.hpp"
#include "algorithms.hpp"
#include "base.hpp"
#include "linear/LaplaceLinear.hpp"
#include "DataVector.h"

namespace sg
{
// ***** OperationBLinear *****

/**
 * Multiplication with vector, not transposed, linear sparse grid
 *
 * @param alpha coefficients of the sparse grid's base functions
 * @param data the vector that should be multiplied
 * @param result the result vector of the matrix vector multiplication
 */
void OperationBLinear::mult(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmB<SLinearBase> op;
	linear_base<unsigned int, unsigned int> base;

	op.mult(storage, base, alpha, data, result);
}

/**
 * Multiplication with vector, transposed, linear sparse grid
 *
 * @param alpha coefficients of the sparse grid's base functions
 * @param data the vector that should be multiplied
 * @param result the result vector of the matrix vector multiplication
 */
void OperationBLinear::multTranspose(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmB<SLinearBase> op;
	linear_base<unsigned int, unsigned int> base;

	op.mult_transpose(storage, base, alpha, data, result);
}

// ***** OperationEvalLinear *****

double OperationEvalLinear::eval(DataVector& alpha, std::vector<double>& point)
{
	typedef std::vector<std::pair<size_t, double> > IndexValVector;

	IndexValVector vec;
	linear_base<unsigned int, unsigned int> base;
	GetAffectedBasisFunctions<linear_base<unsigned int, unsigned int> > ga(storage);

	ga(base, point, vec);

	double result = 0.0;

	for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
	{
		result += iter->second * alpha[iter->first];
	}

	return result;
}

double OperationEvalLinear::test(DataVector& alpha, DataVector& data, DataVector& classes)
{
	linear_base<unsigned int, unsigned int> base;
	return test_dataset(this->storage, base, alpha, data, classes);
}

}
