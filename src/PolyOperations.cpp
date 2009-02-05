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

// Include all operations on the poly grid
#include "basis/poly/operation/OperationBPoly.hpp"
#include "basis/poly/operation/OperationEvalPoly.hpp"
#include "basis/poly/operation/OperationHierarchisationPoly.hpp"

#include "sgpp.hpp"
#include "algorithms.hpp"
#include "basis/basis.hpp"

#include "data/DataVector.h"

#include "exception/operation_exception.hpp"

namespace sg
{

// ***** OperationBPoly *****

void OperationBPoly::mult(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmB<SPolyBase> op;

	op.mult(storage, base, alpha, data, result);
}

void OperationBPoly::multTranspose(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmB<SPolyBase> op;

	op.mult_transpose(storage, base, alpha, data, result);
}

// ***** OperationEvalPoly *****

double OperationEvalPoly::eval(DataVector& alpha, std::vector<double>& point)
{
	typedef std::vector<std::pair<size_t, double> > IndexValVector;

	IndexValVector vec;
	GetAffectedBasisFunctions<SPolyBase> ga(storage);

	ga(base, point, vec);

	double result = 0.0;

	for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
	{
		result += iter->second * alpha[iter->first];
	}

	return result;
}

double OperationEvalPoly::test(DataVector& alpha, DataVector& data, DataVector& classes)
{

	return test_dataset(this->storage, base, alpha, data, classes);
}

// ***** OperationhierarchisationPoly *****

/**
 * Implements the hierarchisation on a sprase grid with poly base functions
 *
 * @param node_values the functions values in the node base
 *
 * @todo Implement the hierarchisation on the sparse grid with poly base functions
 */
void OperationHierarchisationPoly::doHierarchisation(DataVector& node_values)
{
	throw new operation_exception("This operation is not implemented, yet! Sorry ;-)");
}

/**
 * Implements the dehierarchisation on a sprase grid with poly base functions
 *
 * @param alpha the coefficients of the sparse grid's base functions
 *
 * @todo Implement the dehierarchisation on the sparse grid with poly base functions
 */
void OperationHierarchisationPoly::doDehierarchisation(DataVector& alpha)
{
	throw new operation_exception("This operation is not implemented, yet! Sorry ;-)");
}

}


