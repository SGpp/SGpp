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

#include "modlinear/ModLinearOperations.hpp"

#include "sgpp.hpp"
#include "algorithms.hpp"
#include "base.hpp"

#include "modlinear/LaplaceModLinear.hpp"

#include "DataVector.h"

namespace sg
{
// ***** OperationBModLinear *****

/**
 * Multiplication with vector, not transposed, modified linear sparse grid
 *
 * @param alpha coefficients of the sparse grid's base functions
 * @param data the vector that should be multiplied
 * @param result the result vector of the matrix vector multiplication
 */
void OperationBModLinear::mult(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmB<SModLinearBase> op;
	modified_linear_base<unsigned int, unsigned int> base;

	op.mult(storage, base, alpha, data, result);
}

/**
 * Multiplication with vector, transposed, modified linear sparse grid
 *
 * @param alpha coefficients of the sparse grid's base functions
 * @param data the vector that should be multiplied
 * @param result the result vector of the matrix vector multiplication
 */
void OperationBModLinear::multTranspose(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmB<SModLinearBase> op;
	modified_linear_base<unsigned int, unsigned int> base;

	op.mult_transpose(storage, base, alpha, data, result);
}

// ***** OperationEvalModLinear *****

double OperationEvalModLinear::eval(DataVector& alpha, std::vector<double>& point)
{
	typedef std::vector<std::pair<size_t, double> > IndexValVector;

	IndexValVector vec;
	modified_linear_base<unsigned int, unsigned int> base;
	GetAffectedBasisFunctions<modified_linear_base<unsigned int, unsigned int> > ga(storage);

	ga(base, point, vec);

	double result = 0.0;

	for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
	{
		result += iter->second * alpha[iter->first];
	}

	return result;
}

double OperationEvalModLinear::test(DataVector& alpha, DataVector& data, DataVector& classes)
{
	modified_linear_base<unsigned int, unsigned int> base;
	return test_dataset(this->storage, base, alpha, data, classes);
}

// ***** OperationhierarchisationModLinear *****

/**
 * Implements the hierarchisation on a sprase grid with mod linear base functions
 *
 * @param alpha the coefficients of the sparse grid's base functions
 * @param node_values the functions values in the node base
 * @param bDirection this parameter specifies the direction of the hierarchisation: true = hierarchisation, false = dehierarchisation
 *
 * @todo Implement the hierarchisation on the sparse grid with mod linear base functions
 */
void OperationHierarchisationModLinear::PrivateHierarchisation(DataVector& alpha, DataVector& node_values, bool bDirection)
{
	return;
}

/**
 * Implements the hierarchisation on a sprase grid with mod linear base functions
 *
 * @param alpha the coefficients of the sparse grid's base functions
 * @param node_values the functions values in the node base
 *
 * @todo Implement the hierarchisation on the sparse grid with mod linear base functions
 */
void OperationHierarchisationModLinear::doHierarchisation(DataVector& alpha, DataVector& node_values)
{
	PrivateHierarchisation(alpha, node_values, true);
}

/**
 * Implements the dehierarchisation on a sprase grid with mod linear base functions
 *
 * @param alpha the coefficients of the sparse grid's base functions
 * @param node_values the functions values in the node base
 *
 * @todo Implement the dehierarchisation on the sparse grid with mod linear base functions
 */
void OperationHierarchisationModLinear::doDehierarchisation(DataVector& alpha, DataVector& node_values)
{
	PrivateHierarchisation(alpha, node_values, false);
}

};


