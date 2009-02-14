/*
This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems

Copyright (C) 2008 Deepak Pandey

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

// Include all operations on the mod wavelet gird
#include "basis/modwavelet/operation/OperationBModWavelet.hpp"
#include "basis/modwavelet/operation/OperationEvalModWavelet.hpp"
#include "basis/modwavelet/operation/OperationHierarchisationModWavelet.hpp"

#include "sgpp.hpp"

#include "basis/basis.hpp"

#include "data/DataVector.h"

#include "exception/operation_exception.hpp"

namespace sg
{

// ***** OperationBModWavelet *****

void OperationBModWavelet::mult(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmB<SModWaveletBase> op;
	mod_Wavelet_base<unsigned int, unsigned int> base;

	op.mult(storage, base, alpha, data, result);
}

void OperationBModWavelet::multTranspose(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmB<SModWaveletBase> op;
	mod_Wavelet_base<unsigned int, unsigned int> base;

	op.mult_transpose(storage, base, alpha, data, result);
}

// ***** OperationEvalModWavelet *****

double OperationEvalModWavelet::eval(DataVector& alpha, std::vector<double>& point)
{
	typedef std::vector<std::pair<size_t, double> > IndexValVector;

	IndexValVector vec;
	mod_Wavelet_base<unsigned int, unsigned int> base;
	GetAffectedBasisFunctions<mod_Wavelet_base<unsigned int, unsigned int> > ga(storage);

	ga(base, point, vec);

	double result = 0.0;

	for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
	{
		result += iter->second * alpha[iter->first];
	}

	return result;
}

double OperationEvalModWavelet::test(DataVector& alpha, DataVector& data, DataVector& classes)
{
	mod_Wavelet_base<unsigned int, unsigned int> base;
	return test_dataset(this->storage, base, alpha, data, classes);
}

// ***** OperationhierarchisationModWavelet *****

/**
 * Implements the hierarchisation on a sprase grid with mod wavelets base functions
 *
 * @param node_values the functions values in the node base
 *
 * @todo Implement the hierarchisation on the sparse grid with mod wavelets base functions
 */
void OperationhierarchisationModWavelet::doHierarchisation(DataVector& node_values)
{
	throw new operation_exception("This operation is not implemented, yet! Sorry ;-)");
}

/**
 * Implements the dehierarchisation on a sprase grid with mod wavelets base functions
 *
 * @param alpha the coefficients of the sparse grid's base functions
 *
 * @todo Implement the dehierarchisation on the sparse grid with mod wavelets base functions
 */
void OperationhierarchisationModWavelet::doDehierarchisation(DataVector& alpha)
{
	throw new operation_exception("This operation is not implemented, yet! Sorry ;-)");
}

};


