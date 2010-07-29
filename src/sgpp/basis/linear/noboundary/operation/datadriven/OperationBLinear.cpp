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
#include "data/DataMatrix.hpp"
//#include "common/AlignedMemory.hpp"

#define CHUNKDATAPOINTS 512
#define CHUNKGRIDPOINTS 256

namespace sg
{

void OperationBLinear::mult(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	AlgorithmMultipleEvaluation<SLinearBase> op;
	linear_base<unsigned int, unsigned int> base;

	op.mult(storage, base, alpha, data, result);
}

void OperationBLinear::multTranspose(DataVector& alpha, DataMatrix& data, DataVector& result)
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
    	size_t chunksize = (result_size/omp_get_num_threads())+1;
    	size_t start = chunksize*omp_get_thread_num();
    	size_t end = std::min(start+chunksize, result_size);
#else
    	size_t start = 0;
    	size_t end = result_size;
#endif
		for(size_t c = start; c < end; c+=std::min((size_t)CHUNKDATAPOINTS, (end-c)))
		{
			size_t data_end = std::min((size_t)CHUNKDATAPOINTS+c, end);
#ifdef __ICC
			#pragma ivdep
			#pragma vector aligned
#endif
			for (size_t i = c; i < data_end; i++)
			{
				ptrResult[i] = 0.0;
			}

			for (size_t m = 0; m < storageSize; m+=std::min((size_t)CHUNKGRIDPOINTS, (storageSize-m)))
			{
				size_t grid_end = std::min((size_t)CHUNKGRIDPOINTS+m, storageSize);

				for (size_t i = c; i < data_end; i++)
				{
					for (size_t j = m; j < grid_end; j++)
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
		}
#ifdef USEOMP
	}
#endif
}

}
