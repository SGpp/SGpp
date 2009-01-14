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

#include "modpoly/ModPolyOperations.hpp"

#include "sgpp.hpp"
#include "algorithms.hpp"
#include "base.hpp"

#include "DataVector.h"

namespace sg
{


void OperationBModPoly::mult(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmB<SModPolyBase> op;
		
	op.mult(storage, base, alpha, data, result);
}

void OperationBModPoly::multTranspose(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmB<SModPolyBase> op;
		
	op.mult_transpose(storage, base, alpha, data, result);
}


double OperationEvalModPoly::eval(DataVector& alpha, std::vector<double>& point)
{
	typedef std::vector<std::pair<size_t, double> > IndexValVector;
	
	IndexValVector vec;
	GetAffectedBasisFunctions<SModPolyBase> ga(storage);
	
	ga(base, point, vec);

	double result = 0.0;
	
	for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
	{
		result += iter->second * alpha[iter->first];
	}
	
	return result;
}

double OperationEvalModPoly::test(DataVector& alpha, DataVector& data, DataVector& classes)
{
	
	return test_dataset(this->storage, base, alpha, data, classes);
}




}


