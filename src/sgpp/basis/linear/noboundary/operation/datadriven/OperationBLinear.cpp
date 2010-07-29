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
//#include "common/AlignedMemory.hpp"

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

//    double* level = new double[storage->size()*storage->dim()];
//    double* index = new double[storage->size()*storage->dim()];
//
//    storage->getLevelIndexArraysForEval(level, index);
//    multTransposeIterative(level, index, alpha, data, result);
//
//    delete[] level;
//    delete[] index;
}

void OperationBLinear::multTransposeIterative(double* Level, double* Index, DataVector& alpha, DataVector& data, DataVector& result)
{
	size_t result_size = result.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();
    double* ptrAlpha = alpha.getPointer();
    double* ptrData = data.getPointer();
    double* ptrResult = result.getPointer();

#ifdef USEOMP
	#pragma omp parallel
	{
		#pragma omp for schedule (static)
		for(size_t i = 0; i < result_size; i++)
		{
	        ptrResult[i] = 0.0;

#ifdef __ICC
			#pragma ivdep
			#pragma vector aligned
#endif
	        for (size_t j = 0; j < storageSize; j++)
	        {
	        	double curSupport = 1.0;
#ifdef __ICC
				#pragma ivdep
				#pragma vector aligned
#endif
	        	for (size_t d = 0; d < dims; d++)
	        	{
	        		double eval = (Level[(j*dims)+d]) * (ptrData[(i*dims)+d]);
	        		double index_calc = eval - (Index[(j*dims)+d]);
	        		double abs = fabs(index_calc);
	        		double last = 1.0 - abs;
	        		double localSupport = std::max(last, 0.0);
	        		curSupport *= localSupport;
	        	}

	        	ptrResult[i] += (curSupport * ptrAlpha[j]);
	        }
		}
	}
#else
	for(size_t i = 0; i < result_size; i++)
	{
        ptrResult[i] = 0.0;

#ifdef __ICC
		#pragma ivdep
		#pragma vector aligned
#endif
        for (size_t j = 0; j < storageSize; j++)
        {
        	double curSupport = 1.0;

#ifdef __ICC
			#pragma ivdep
			#pragma vector aligned
#endif
        	for (size_t d = 0; d < dims; d++)
        	{
        		double eval = (Level[(j*dims)+d]) * (ptrData[(i*dims)+d]);
        		double index_calc = eval - (Index[(j*dims)+d]);
        		double abs = fabs(index_calc);
        		double last = 1.0 - abs;
        		double localSupport = std::max(last, 0.0);
        		curSupport *= localSupport;
        	}

        	ptrResult[i] += (curSupport * ptrAlpha[j]);
        }
	}
#endif
}

}
