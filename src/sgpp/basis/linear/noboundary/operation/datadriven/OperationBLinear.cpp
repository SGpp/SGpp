/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"

#include "basis/basis.hpp"

#include "basis/linear/noboundary/operation/datadriven/OperationBLinear.hpp"

#include "data/DataVector.hpp"

namespace sg
{

void OperationBLinear::mult(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmMultipleEvaluation<SLinearBase> op;
	linear_base<unsigned int, unsigned int> base;

	op.mult(storage, base, alpha, data, result);
}

void OperationBLinear::multTranspose(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmMultipleEvaluation<SLinearBase> op;
	linear_base<unsigned int, unsigned int> base;

	op.mult_transpose(storage, base, alpha, data, result);
	//op.mult_transpose_iterative(storage, base, alpha, data, result);
}

void OperationBLinear::multTransposeIterativeTest(unsigned int* Level, unsigned int* Index, DataVector& alpha, DataVector& data, DataVector& result)
{
	size_t result_size = result.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

#ifdef USEOMP
	#pragma omp parallel
	{
#ifdef __ICC
		//#pragma ivdep
#endif
		#pragma omp for schedule (static)
		for(size_t i = 0; i < result_size; i++)
		{
	        result[i] = 0.0;
	        double curSupport;

#ifdef __ICC
			//#pragma ivdep
			//#pragma vector aligend
#endif
	        for (size_t j = 0; j < storageSize; j++)
	        {
	        	curSupport = 1.0;

#ifdef __ICC
				//#pragma ivdep
				//#pragma vector aligend
#endif
	        	for (size_t d = 0; d < dims; d++)
	        	{
	        		curSupport *= std::max(1.0 - fabs(((1<<(Level[(j*dims)+d])) * (data[(i*dims)+d])) - (Index[(j*dims)+d])), 0.0);
	        	}

	        	result[i] += (curSupport * alpha[j]);
	        }
		}
	}
#else
	for(size_t i = 0; i < result_size; i++)
	{
        result[i] = 0.0;
        double curSupport;

        for (size_t j = 0; j < storageSize; j++)
        {
        	curSupport = 1.0;

        	for (size_t d = 0; d < dims; d++)
        	{
        		curSupport *= std::max(1.0 - fabs(((1<<(Level[(j*dims)+d])) * (data[(i*dims)+d])) - (Index[(j*dims)+d])), 0.0);
        	}

        	result[i] += (curSupport * alpha[j]);
        }
	}
#endif
}

}
