/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeSSELinear.hpp"
#include "exception/operation_exception.hpp"

#ifdef USEOMP
#include "omp.h"
#endif

#ifdef USEICCINTRINSICS
// include SSE3 intrinsics
#include <pmmintrin.h>

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

#define CHUNKDATAPOINTS 12 // must be divide-able by 12
#define CHUNKGRIDPOINTS 12

namespace sg
{

OperationBIterativeSSELinear::OperationBIterativeSSELinear(GridStorage* storage) : storage(storage)
{
	Level = new DataMatrix(storage->size(), storage->dim());
	Index = new DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*Level, *Index);
}

OperationBIterativeSSELinear::~OperationBIterativeSSELinear()
{
	delete Level;
	delete Index;
}

void OperationBIterativeSSELinear::rebuildLevelAndIndex()
{
	delete Level;
	delete Index;

	Level = new DataMatrix(storage->size(), storage->dim());
	Index = new DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*Level, *Index);
}

void OperationBIterativeSSELinear::multVectorized(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	size_t source_size = alpha.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();
    double* ptrSource = alpha.getPointer();
    double* ptrData = data.getPointer();
    double* ptrLevel = this->Level->getPointer();
    double* ptrIndex = this->Index->getPointer();

    if (data.getNcols() % 2 != 0 || source_size != data.getNcols())
    {
    	throw operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    	return;
    }

    result.setAll(0.0);

#ifdef USEOMP
    #pragma omp parallel
	{
		size_t chunksize = (source_size/omp_get_num_threads())+1;
		// assure that every subarray is 16-byte aligned
		if (chunksize % 2 != 0)
		{
			chunksize++;
		}
    	size_t start = chunksize*omp_get_thread_num();
    	size_t end = std::min<size_t>(start+chunksize, source_size);

    	DataVector myResult(result.getSize());
    	myResult.setAll(0.0);
    	double* ptrResult = myResult.getPointer();
#else
    	size_t start = 0;
    	size_t end = source_size;
    	double* ptrResult = result.getPointer();
#endif
		for(size_t c = start; c < end; c+=std::min<size_t>((size_t)CHUNKDATAPOINTS, (end-c)))
		{
			size_t data_end = std::min<size_t>((size_t)CHUNKDATAPOINTS+c, end);

			for (size_t m = 0; m < storageSize; m+=std::min<size_t>((size_t)CHUNKGRIDPOINTS, (storageSize-m)))
			{
				size_t grid_end = std::min<size_t>((size_t)CHUNKGRIDPOINTS+m, storageSize);
#ifdef USEICCINTRINSICS
				if ((data_end-c) == CHUNKDATAPOINTS && (grid_end-m) == CHUNKGRIDPOINTS)
				{
					for (size_t i = c; i < c+CHUNKDATAPOINTS; i+=12)
					{
						for (size_t j = m; j < m+CHUNKGRIDPOINTS; j++)
						{
							__m128d support_0 = _mm_load_pd(&(ptrSource[i]));
							__m128d support_1 = _mm_load_pd(&(ptrSource[i+2]));
							__m128d support_2 = _mm_load_pd(&(ptrSource[i+4]));
							__m128d support_3 = _mm_load_pd(&(ptrSource[i+6]));
							__m128d support_4 = _mm_load_pd(&(ptrSource[i+8]));
							__m128d support_5 = _mm_load_pd(&(ptrSource[i+10]));

							__m128d one = _mm_set1_pd(1.0);
							__m128d zero = _mm_set1_pd(0.0);

							for (size_t d = 0; d < dims; d++)
							{
								__m128d eval_0 = _mm_load_pd(&(ptrData[(d*source_size)+i]));
								__m128d eval_1 = _mm_load_pd(&(ptrData[(d*source_size)+i+2]));
								__m128d eval_2 = _mm_load_pd(&(ptrData[(d*source_size)+i+4]));
								__m128d eval_3 = _mm_load_pd(&(ptrData[(d*source_size)+i+6]));
								__m128d eval_4 = _mm_load_pd(&(ptrData[(d*source_size)+i+8]));
								__m128d eval_5 = _mm_load_pd(&(ptrData[(d*source_size)+i+10]));

								__m128d level = _mm_loaddup_pd(&(ptrLevel[(j*dims)+d]));
								__m128d index = _mm_loaddup_pd(&(ptrIndex[(j*dims)+d]));

								eval_0 = _mm_mul_pd(eval_0, level);
								eval_1 = _mm_mul_pd(eval_1, level);
								eval_2 = _mm_mul_pd(eval_2, level);
								eval_3 = _mm_mul_pd(eval_3, level);
								eval_4 = _mm_mul_pd(eval_4, level);
								eval_5 = _mm_mul_pd(eval_5, level);

								eval_0 = _mm_sub_pd(eval_0, index);
								eval_1 = _mm_sub_pd(eval_1, index);
								eval_2 = _mm_sub_pd(eval_2, index);
								eval_3 = _mm_sub_pd(eval_3, index);
								eval_4 = _mm_sub_pd(eval_4, index);
								eval_5 = _mm_sub_pd(eval_5, index);

								eval_0 = _mm_abs_pd(eval_0);
								eval_1 = _mm_abs_pd(eval_1);
								eval_2 = _mm_abs_pd(eval_2);
								eval_3 = _mm_abs_pd(eval_3);
								eval_4 = _mm_abs_pd(eval_4);
								eval_5 = _mm_abs_pd(eval_5);

								eval_0 = _mm_sub_pd(one, eval_0);
								eval_1 = _mm_sub_pd(one, eval_1);
								eval_2 = _mm_sub_pd(one, eval_2);
								eval_3 = _mm_sub_pd(one, eval_3);
								eval_4 = _mm_sub_pd(one, eval_4);
								eval_5 = _mm_sub_pd(one, eval_5);

								eval_0 = _mm_max_pd(zero, eval_0);
								eval_1 = _mm_max_pd(zero, eval_1);
								eval_2 = _mm_max_pd(zero, eval_2);
								eval_3 = _mm_max_pd(zero, eval_3);
								eval_4 = _mm_max_pd(zero, eval_4);
								eval_5 = _mm_max_pd(zero, eval_5);

								support_0 = _mm_mul_pd(support_0, eval_0);
								support_1 = _mm_mul_pd(support_1, eval_1);
								support_2 = _mm_mul_pd(support_2, eval_2);
								support_3 = _mm_mul_pd(support_3, eval_3);
								support_4 = _mm_mul_pd(support_4, eval_4);
								support_5 = _mm_mul_pd(support_5, eval_5);
							}

							__m128d res_0 = _mm_setzero_pd();
							res_0 = _mm_loadl_pd(res_0, &(ptrResult[j]));

							support_0 = _mm_add_pd(support_0, support_1);
							support_2 = _mm_add_pd(support_2, support_3);
							support_4 = _mm_add_pd(support_4, support_5);
							support_0 = _mm_add_pd(support_0, support_2);
							support_0 = _mm_add_pd(support_0, support_4);

							support_0 = _mm_hadd_pd(support_0, support_0);
							res_0 = _mm_add_sd(res_0, support_0);

							_mm_storel_pd(&(ptrResult[j]), res_0);
						}
					}
				}
				else
				{
					for (size_t i = c; i < data_end; i++)
					{
						for (size_t j = m; j < grid_end; j++)
						{
							double curSupport = ptrSource[i];

							#pragma ivdep
							#pragma vector aligned
							for (size_t d = 0; d < dims; d++)
							{
								double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*source_size)+i]));
								double index_calc = eval - (ptrIndex[(j*dims)+d]);
								double abs = fabs(index_calc);
								double last = 1.0 - abs;
								double localSupport = std::max<double>(last, 0.0);
								curSupport *= localSupport;
							}

							ptrResult[j] += curSupport;
						}
					}
				}
#else
				for (size_t i = c; i < data_end; i++)
				{
					for (size_t j = m; j < grid_end; j++)
					{
						double curSupport = ptrSource[i];
#ifdef __ICC
						#pragma ivdep
						#pragma vector aligned
#endif
						for (size_t d = 0; d < dims; d++)
						{
							double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*source_size)+i]));
							double index_calc = eval - (ptrIndex[(j*dims)+d]);
							double abs = fabs(index_calc);
							double last = 1.0 - abs;
							double localSupport = std::max<double>(last, 0.0);
							curSupport *= localSupport;
						}

						ptrResult[j] += curSupport;
					}
				}
#endif
	        }
		}
#ifdef USEOMP
		// sum private result vectors
		#pragma omp critical
		{
			result.add(myResult);
		}
	}
#endif
}

void OperationBIterativeSSELinear::multTransposeVectorized(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	size_t result_size = result.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();
    double* ptrAlpha = alpha.getPointer();
    double* ptrData = data.getPointer();
    double* ptrResult = result.getPointer();
    double* ptrLevel = this->Level->getPointer();
    double* ptrIndex = this->Index->getPointer();

    if (data.getNcols() % 2 != 0 || result_size != data.getNcols())
    {
    	throw operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    	return;
    }

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
    	size_t end = std::min<size_t>(start+chunksize, result_size);
#else
    	size_t start = 0;
    	size_t end = result_size;
#endif
		for(size_t c = start; c < end; c+=std::min<size_t>((size_t)CHUNKDATAPOINTS, (end-c)))
		{
			size_t data_end = std::min<size_t>((size_t)CHUNKDATAPOINTS+c, end);

#ifdef __ICC
			#pragma ivdep
			#pragma vector aligned
#endif
			for (size_t i = c; i < data_end; i++)
			{
				ptrResult[i] = 0.0;
			}

			for (size_t m = 0; m < storageSize; m+=std::min<size_t>((size_t)CHUNKGRIDPOINTS, (storageSize-m)))
			{
				size_t grid_end = std::min<size_t>((size_t)CHUNKGRIDPOINTS+m, storageSize);
#ifdef USEICCINTRINSICS
				if ((data_end-c) == CHUNKDATAPOINTS && (grid_end-m) == CHUNKGRIDPOINTS)
				{
					for (size_t i = c; i < c+CHUNKDATAPOINTS; i+=12)
					{
						for (size_t j = m; j < m+CHUNKGRIDPOINTS; j++)
						{
							__m128d support_0 = _mm_loaddup_pd(&(ptrAlpha[j]));
							__m128d support_1 = _mm_loaddup_pd(&(ptrAlpha[j]));
							__m128d support_2 = _mm_loaddup_pd(&(ptrAlpha[j]));
							__m128d support_3 = _mm_loaddup_pd(&(ptrAlpha[j]));
							__m128d support_4 = _mm_loaddup_pd(&(ptrAlpha[j]));
							__m128d support_5 = _mm_loaddup_pd(&(ptrAlpha[j]));

							__m128d one = _mm_set1_pd(1.0);
							__m128d zero = _mm_set1_pd(0.0);

							for (size_t d = 0; d < dims; d++)
							{
								__m128d eval_0 = _mm_load_pd(&(ptrData[(d*result_size)+i]));
								__m128d eval_1 = _mm_load_pd(&(ptrData[(d*result_size)+i+2]));
								__m128d eval_2 = _mm_load_pd(&(ptrData[(d*result_size)+i+4]));
								__m128d eval_3 = _mm_load_pd(&(ptrData[(d*result_size)+i+6]));
								__m128d eval_4 = _mm_load_pd(&(ptrData[(d*result_size)+i+8]));
								__m128d eval_5 = _mm_load_pd(&(ptrData[(d*result_size)+i+10]));

								__m128d level = _mm_loaddup_pd(&(ptrLevel[(j*dims)+d]));
								__m128d index = _mm_loaddup_pd(&(ptrIndex[(j*dims)+d]));

								eval_0 = _mm_mul_pd(eval_0, level);
								eval_1 = _mm_mul_pd(eval_1, level);
								eval_2 = _mm_mul_pd(eval_2, level);
								eval_3 = _mm_mul_pd(eval_3, level);
								eval_4 = _mm_mul_pd(eval_4, level);
								eval_5 = _mm_mul_pd(eval_5, level);

								eval_0 = _mm_sub_pd(eval_0, index);
								eval_1 = _mm_sub_pd(eval_1, index);
								eval_2 = _mm_sub_pd(eval_2, index);
								eval_3 = _mm_sub_pd(eval_3, index);
								eval_4 = _mm_sub_pd(eval_4, index);
								eval_5 = _mm_sub_pd(eval_5, index);

								eval_0 = _mm_abs_pd(eval_0);
								eval_1 = _mm_abs_pd(eval_1);
								eval_2 = _mm_abs_pd(eval_2);
								eval_3 = _mm_abs_pd(eval_3);
								eval_4 = _mm_abs_pd(eval_4);
								eval_5 = _mm_abs_pd(eval_5);

								eval_0 = _mm_sub_pd(one, eval_0);
								eval_1 = _mm_sub_pd(one, eval_1);
								eval_2 = _mm_sub_pd(one, eval_2);
								eval_3 = _mm_sub_pd(one, eval_3);
								eval_4 = _mm_sub_pd(one, eval_4);
								eval_5 = _mm_sub_pd(one, eval_5);

								eval_0 = _mm_max_pd(zero, eval_0);
								eval_1 = _mm_max_pd(zero, eval_1);
								eval_2 = _mm_max_pd(zero, eval_2);
								eval_3 = _mm_max_pd(zero, eval_3);
								eval_4 = _mm_max_pd(zero, eval_4);
								eval_5 = _mm_max_pd(zero, eval_5);

								support_0 = _mm_mul_pd(support_0, eval_0);
								support_1 = _mm_mul_pd(support_1, eval_1);
								support_2 = _mm_mul_pd(support_2, eval_2);
								support_3 = _mm_mul_pd(support_3, eval_3);
								support_4 = _mm_mul_pd(support_4, eval_4);
								support_5 = _mm_mul_pd(support_5, eval_5);
							}

							__m128d res_0 = _mm_load_pd(&(ptrResult[i]));
							__m128d res_1 = _mm_load_pd(&(ptrResult[i+2]));
							__m128d res_2 = _mm_load_pd(&(ptrResult[i+4]));
							__m128d res_3 = _mm_load_pd(&(ptrResult[i+6]));
							__m128d res_4 = _mm_load_pd(&(ptrResult[i+8]));
							__m128d res_5 = _mm_load_pd(&(ptrResult[i+10]));

							res_0 = _mm_add_pd(res_0, support_0);
							res_1 = _mm_add_pd(res_1, support_1);
							res_2 = _mm_add_pd(res_2, support_2);
							res_3 = _mm_add_pd(res_3, support_3);
							res_4 = _mm_add_pd(res_4, support_4);
							res_5 = _mm_add_pd(res_5, support_5);

							_mm_store_pd(&(ptrResult[i]), res_0);
							_mm_store_pd(&(ptrResult[i+2]), res_1);
							_mm_store_pd(&(ptrResult[i+4]), res_2);
							_mm_store_pd(&(ptrResult[i+6]), res_3);
							_mm_store_pd(&(ptrResult[i+8]), res_4);
							_mm_store_pd(&(ptrResult[i+10]), res_5);
						}
					}
				}
				else
				{
					for (size_t i = c; i < data_end; i++)
					{
						for (size_t j = m; j < grid_end; j++)
						{
							double curSupport = ptrAlpha[j];

							#pragma ivdep
							#pragma vector aligned
							for (size_t d = 0; d < dims; d++)
							{
								double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*result_size)+i]));
								double index_calc = eval - (ptrIndex[(j*dims)+d]);
								double abs = fabs(index_calc);
								double last = 1.0 - abs;
								double localSupport = std::max<double>(last, 0.0);
								curSupport *= localSupport;
							}

							ptrResult[i] += curSupport;
						}
					}
				}
#else
				for (size_t i = c; i < data_end; i++)
				{
					for (size_t j = m; j < grid_end; j++)
					{
						double curSupport = ptrAlpha[j];
#ifdef __ICC
						#pragma ivdep
						#pragma vector aligned
#endif
						for (size_t d = 0; d < dims; d++)
						{
							double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*result_size)+i]));
							double index_calc = eval - (ptrIndex[(j*dims)+d]);
							double abs = fabs(index_calc);
							double last = 1.0 - abs;
							double localSupport = std::max<double>(last, 0.0);
							curSupport *= localSupport;
						}

						ptrResult[i] += curSupport;
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
