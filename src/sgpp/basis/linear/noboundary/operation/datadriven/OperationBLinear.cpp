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

#ifdef USEINTRINSICS
#include <xmmintrin.h>

union doubleAbsMask
{
   const double d;
   const __int64 i;

   doubleAbsMask() : i(0x7FFFFFFFFFFFFFFF) {}
};

_MM_ALIGN16 const doubleAbsMask absMask;
static const __m128d abs2Mask = _mm_load1_pd( &absMask.d );

const __m128d _mm_abs_pd( const __m128d& x)
{
       return _mm_and_pd( abs2Mask, x);
}
#endif

#define CHUNKDATAPOINTS 256 // must be divide-able by 8
#define CHUNKGRIDPOINTS 128

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
		size_t chunksize = (result_size/omp_get_num_threads())+1;
		// assure that every subarray is 16-byte aligned
		if (chunksize % 2 != 0)
		{
			chunksize++;
		}
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

#ifdef USEINTRINSICS
				if ((data_end-c) == CHUNKDATAPOINTS && (grid_end-m) == CHUNKGRIDPOINTS)
				{
					for (size_t i = c; i < c+CHUNKDATAPOINTS; i+=8)
					{
						for (size_t j = m; j < m+CHUNKGRIDPOINTS; j++)
						{
							__m128d support_0 = _mm_set1_pd(1.0);
							__m128d support_1 = _mm_set1_pd(1.0);
							__m128d support_2 = _mm_set1_pd(1.0);
							__m128d support_3 = _mm_set1_pd(1.0);

							for (size_t d = 0; d < dims; d++)
							{
								// do data transpose before
								__m128d eval_0 = _mm_set_pd(ptrData[((i+1)*dims)+d], ptrData[(i*dims)+d]);
								__m128d eval_1 = _mm_set_pd(ptrData[((i+3)*dims)+d], ptrData[((i+2)*dims)+d]);
								__m128d eval_2 = _mm_set_pd(ptrData[((i+5)*dims)+d], ptrData[((i+4)*dims)+d]);
								__m128d eval_3 = _mm_set_pd(ptrData[((i+7)*dims)+d], ptrData[((i+6)*dims)+d]);

								__m128d level = _mm_loaddup_pd(&(Level[(j*dims)+d]));
								__m128d index = _mm_loaddup_pd(&(Index[(j*dims)+d]));

								eval_0 = _mm_mul_pd(eval_0, level);
								eval_1 = _mm_mul_pd(eval_1, level);
								eval_2 = _mm_mul_pd(eval_2, level);
								eval_3 = _mm_mul_pd(eval_3, level);

								eval_0 = _mm_sub_pd(eval_0, index);
								eval_1 = _mm_sub_pd(eval_1, index);
								eval_2 = _mm_sub_pd(eval_2, index);
								eval_3 = _mm_sub_pd(eval_3, index);

								__m128d one = _mm_set1_pd(1.0);
								__m128d zero = _mm_set1_pd(0.0);

								eval_0 = _mm_abs_pd(eval_0);
								eval_1 = _mm_abs_pd(eval_1);
								eval_2 = _mm_abs_pd(eval_2);
								eval_3 = _mm_abs_pd(eval_3);

								eval_0 = _mm_sub_pd(one, eval_0);
								eval_1 = _mm_sub_pd(one, eval_1);
								eval_2 = _mm_sub_pd(one, eval_2);
								eval_3 = _mm_sub_pd(one, eval_3);

								eval_0 = _mm_max_pd(zero, eval_0);
								eval_1 = _mm_max_pd(zero, eval_1);
								eval_2 = _mm_max_pd(zero, eval_2);
								eval_3 = _mm_max_pd(zero, eval_3);

								support_0 = _mm_mul_pd(support_0, eval_0);
								support_1 = _mm_mul_pd(support_1, eval_1);
								support_2 = _mm_mul_pd(support_2, eval_2);
								support_3 = _mm_mul_pd(support_3, eval_3);
							}

							__m128d alpha = _mm_loaddup_pd(&(ptrAlpha[j]));
							__m128d res_0 = _mm_load_pd(&(ptrResult[i]));
							__m128d res_1 = _mm_load_pd(&(ptrResult[i+2]));
							__m128d res_2 = _mm_load_pd(&(ptrResult[i+4]));
							__m128d res_3 = _mm_load_pd(&(ptrResult[i+6]));

							support_0 = _mm_mul_pd(support_0, alpha);
							support_1 = _mm_mul_pd(support_1, alpha);
							support_2 = _mm_mul_pd(support_2, alpha);
							support_3 = _mm_mul_pd(support_3, alpha);

							res_0 = _mm_add_pd(res_0, support_0);
							res_1 = _mm_add_pd(res_1, support_1);
							res_2 = _mm_add_pd(res_2, support_2);
							res_3 = _mm_add_pd(res_3, support_3);

							_mm_store_pd(&(ptrResult[i]), res_0);
							_mm_store_pd(&(ptrResult[i+2]), res_1);
							_mm_store_pd(&(ptrResult[i+4]), res_2);
							_mm_store_pd(&(ptrResult[i+6]), res_3);
						}
					}
				}
				else
				{
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
								double eval = ((Level[(j*dims)+d]) * (ptrData[(i*dims)+d]));
								double index_calc = eval - (Index[(j*dims)+d]);
								double abs = fabs(index_calc);
								double last = 1.0 - abs;
								double localSupport = std::max(last, 0.0);
								curSupport *= localSupport;

	//							curSupport *= std::max(1.0 - fabs(((Level[(j*dims)+d]) * (ptrData[(i*dims)+d])) - (Index[(j*dims)+d])), 0.0);
							}

							ptrResult[i] += (curSupport * ptrAlpha[j]);
						}
					}
				}
#else
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
						double eval = ((Level[(j*dims)+d]) * (ptrData[(i*dims)+d]));
						double index_calc = eval - (Index[(j*dims)+d]);
						double abs = fabs(index_calc);
						double last = 1.0 - abs;
						double localSupport = std::max(last, 0.0);
						curSupport *= localSupport;

//							curSupport *= std::max(1.0 - fabs(((Level[(j*dims)+d]) * (ptrData[(i*dims)+d])) - (Index[(j*dims)+d])), 0.0);
					}

					ptrResult[i] += (curSupport * ptrAlpha[j]);
				}
			}
#endif

	        }
		}
#ifdef USEOMP
	}
#endif
}

}
