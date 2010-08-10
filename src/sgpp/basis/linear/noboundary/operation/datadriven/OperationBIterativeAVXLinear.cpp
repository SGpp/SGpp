/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeAVXLinear.hpp"
#include "exception/operation_exception.hpp"

#ifdef USEOMP
#include "omp.h"
#endif

#ifdef USEICCINTRINSICS
#ifdef USEAVX
#include <immintrin.h>
#else
#include "common/avxintrin_emu.h"
#endif

union doubleAbsMaskAVX
{
   const double d;
   const __int64 i;

   doubleAbsMaskAVX() : i(0x7FFFFFFFFFFFFFFF) {}
};

_MM_ALIGN16 const doubleAbsMaskAVX absMaskAVX;

static const __m256d abs2MaskAVX = _mm256_broadcast_sd( &(absMaskAVX.d) );

const __m256d _mm256_abs_pd( const __m256d& x)
{
       return _mm256_and_pd( abs2MaskAVX, x);
}

static const __m256i ldStMaskAVX = _mm256_set_epi64x(0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xFFFFFFFFFFFFFFFF);
#endif

#define CHUNKDATAPOINTS_AVX 24 // must be divide-able by 24
#define CHUNKGRIDPOINTS_AVX 12

namespace sg
{

OperationBIterativeAVXLinear::OperationBIterativeAVXLinear(GridStorage* storage) : storage(storage)
{
	Level = new DataMatrix(storage->size(), storage->dim());
	Index = new DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*Level, *Index);
}

OperationBIterativeAVXLinear::~OperationBIterativeAVXLinear()
{
	delete Level;
	delete Index;
}

void OperationBIterativeAVXLinear::rebuildLevelAndIndex()
{
	delete Level;
	delete Index;

	Level = new DataMatrix(storage->size(), storage->dim());
	Index = new DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*Level, *Index);
}

void OperationBIterativeAVXLinear::multVectorized(DataVector& alpha, DataMatrix& data, DataVector& result)
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
    	size_t end = std::min(start+chunksize, source_size);

    	DataVector myResult(result.getSize());
    	myResult.setAll(0.0);
    	double* ptrResult = myResult.getPointer();
#else
    	size_t start = 0;
    	size_t end = source_size;
    	double* ptrResult = result.getPointer();
#endif
		for(size_t c = start; c < end; c+=std::min((size_t)CHUNKDATAPOINTS_AVX, (end-c)))
		{
			size_t data_end = std::min((size_t)CHUNKDATAPOINTS_AVX+c, end);

			for (size_t m = 0; m < storageSize; m+=std::min((size_t)CHUNKGRIDPOINTS_AVX, (storageSize-m)))
			{
				size_t grid_end = std::min((size_t)CHUNKGRIDPOINTS_AVX+m, storageSize);
#ifdef USEICCINTRINSICS
				if ((data_end-c) == CHUNKDATAPOINTS_AVX && (grid_end-m) == CHUNKGRIDPOINTS_AVX)
				{
					for (size_t i = c; i < c+CHUNKDATAPOINTS_AVX; i+=24)
					{
						for (size_t j = m; j < m+CHUNKGRIDPOINTS_AVX; j++)
						{
							__m256d support_0 = _mm256_load_pd(&(ptrSource[i]));
							__m256d support_1 = _mm256_load_pd(&(ptrSource[i+4]));
							__m256d support_2 = _mm256_load_pd(&(ptrSource[i+8]));
							__m256d support_3 = _mm256_load_pd(&(ptrSource[i+12]));
							__m256d support_4 = _mm256_load_pd(&(ptrSource[i+16]));
							__m256d support_5 = _mm256_load_pd(&(ptrSource[i+20]));

							__m256d one = _mm256_set1_pd(1.0);
							__m256d zero = _mm256_set1_pd(0.0);

							for (size_t d = 0; d < dims; d++)
							{
								__m256d eval_0 = _mm256_load_pd(&(ptrData[(d*source_size)+i]));
								__m256d eval_1 = _mm256_load_pd(&(ptrData[(d*source_size)+i+4]));
								__m256d eval_2 = _mm256_load_pd(&(ptrData[(d*source_size)+i+8]));
								__m256d eval_3 = _mm256_load_pd(&(ptrData[(d*source_size)+i+12]));
								__m256d eval_4 = _mm256_load_pd(&(ptrData[(d*source_size)+i+16]));
								__m256d eval_5 = _mm256_load_pd(&(ptrData[(d*source_size)+i+20]));

								__m256d level = _mm256_broadcast_sd(&(ptrLevel[(j*dims)+d]));
								__m256d index = _mm256_broadcast_sd(&(ptrIndex[(j*dims)+d]));

								eval_0 = _mm256_mul_pd(eval_0, level);
								eval_1 = _mm256_mul_pd(eval_1, level);
								eval_2 = _mm256_mul_pd(eval_2, level);
								eval_3 = _mm256_mul_pd(eval_3, level);
								eval_4 = _mm256_mul_pd(eval_4, level);
								eval_5 = _mm256_mul_pd(eval_5, level);

								eval_0 = _mm256_sub_pd(eval_0, index);
								eval_1 = _mm256_sub_pd(eval_1, index);
								eval_2 = _mm256_sub_pd(eval_2, index);
								eval_3 = _mm256_sub_pd(eval_3, index);
								eval_4 = _mm256_sub_pd(eval_4, index);
								eval_5 = _mm256_sub_pd(eval_5, index);

								eval_0 = _mm256_abs_pd(eval_0);
								eval_1 = _mm256_abs_pd(eval_1);
								eval_2 = _mm256_abs_pd(eval_2);
								eval_3 = _mm256_abs_pd(eval_3);
								eval_4 = _mm256_abs_pd(eval_4);
								eval_5 = _mm256_abs_pd(eval_5);

								eval_0 = _mm256_sub_pd(one, eval_0);
								eval_1 = _mm256_sub_pd(one, eval_1);
								eval_2 = _mm256_sub_pd(one, eval_2);
								eval_3 = _mm256_sub_pd(one, eval_3);
								eval_4 = _mm256_sub_pd(one, eval_4);
								eval_5 = _mm256_sub_pd(one, eval_5);

								eval_0 = _mm256_max_pd(zero, eval_0);
								eval_1 = _mm256_max_pd(zero, eval_1);
								eval_2 = _mm256_max_pd(zero, eval_2);
								eval_3 = _mm256_max_pd(zero, eval_3);
								eval_4 = _mm256_max_pd(zero, eval_4);
								eval_5 = _mm256_max_pd(zero, eval_5);

								support_0 = _mm256_mul_pd(support_0, eval_0);
								support_1 = _mm256_mul_pd(support_1, eval_1);
								support_2 = _mm256_mul_pd(support_2, eval_2);
								support_3 = _mm256_mul_pd(support_3, eval_3);
								support_4 = _mm256_mul_pd(support_4, eval_4);
								support_5 = _mm256_mul_pd(support_5, eval_5);
							}

							__m256d res_0 = _mm256_maskload_pd(&(ptrResult[j]), ldStMaskAVX);

							support_0 = _mm256_add_pd(support_0, support_1);
							support_2 = _mm256_add_pd(support_2, support_3);
							support_4 = _mm256_add_pd(support_4, support_5);
							support_0 = _mm256_add_pd(support_0, support_2);
							support_0 = _mm256_add_pd(support_0, support_4);

							support_0 = _mm256_hadd_pd(support_0, support_0);
							__m256d tmp = _mm256_permute2f128_pd(support_0, support_0, 0x81);
							support_0 = _mm256_add_pd(support_0, tmp);
							res_0 = _mm256_add_pd(res_0, support_0);

							_mm256_maskstore_pd(&(ptrResult[j]), ldStMaskAVX, res_0);
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

							#pragma ivdep
							#pragma vector aligned
							for (size_t d = 0; d < dims; d++)
							{
								double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*source_size)+i]));
								double index_calc = eval - (ptrIndex[(j*dims)+d]);
								double abs = fabs(index_calc);
								double last = 1.0 - abs;
								double localSupport = std::max(last, 0.0);
								curSupport *= localSupport;
							}

							ptrResult[j] += (curSupport * ptrSource[i]);
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
							double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*source_size)+i]));
							double index_calc = eval - (ptrIndex[(j*dims)+d]);
							double abs = fabs(index_calc);
							double last = 1.0 - abs;
							double localSupport = std::max(last, 0.0);
							curSupport *= localSupport;
						}

						ptrResult[j] += (curSupport * ptrSource[i]);
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

void OperationBIterativeAVXLinear::multTransposeVectorized(DataVector& alpha, DataMatrix& data, DataVector& result)
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
    	size_t end = std::min(start+chunksize, result_size);
#else
    	size_t start = 0;
    	size_t end = result_size;
#endif
		for(size_t c = start; c < end; c+=std::min((size_t)CHUNKDATAPOINTS_AVX, (end-c)))
		{
			size_t data_end = std::min((size_t)CHUNKDATAPOINTS_AVX+c, end);

#ifdef __ICC
			#pragma ivdep
			#pragma vector aligned
#endif
			for (size_t i = c; i < data_end; i++)
			{
				ptrResult[i] = 0.0;
			}

			for (size_t m = 0; m < storageSize; m+=std::min((size_t)CHUNKGRIDPOINTS_AVX, (storageSize-m)))
			{
				size_t grid_end = std::min((size_t)CHUNKGRIDPOINTS_AVX+m, storageSize);
#ifdef USEICCINTRINSICS
				if ((data_end-c) == CHUNKDATAPOINTS_AVX && (grid_end-m) == CHUNKGRIDPOINTS_AVX)
				{
					for (size_t i = c; i < c+CHUNKDATAPOINTS_AVX; i+=24)
					{
						for (size_t j = m; j < m+CHUNKGRIDPOINTS_AVX; j++)
						{
							__m256d support_0 = _mm256_broadcast_sd(&(ptrAlpha[j]));
							__m256d support_1 = _mm256_broadcast_sd(&(ptrAlpha[j]));
							__m256d support_2 = _mm256_broadcast_sd(&(ptrAlpha[j]));
							__m256d support_3 = _mm256_broadcast_sd(&(ptrAlpha[j]));
							__m256d support_4 = _mm256_broadcast_sd(&(ptrAlpha[j]));
							__m256d support_5 = _mm256_broadcast_sd(&(ptrAlpha[j]));

							__m256d one = _mm256_set1_pd(1.0);
							__m256d zero = _mm256_set1_pd(0.0);

							for (size_t d = 0; d < dims; d++)
							{
								__m256d eval_0 = _mm256_load_pd(&(ptrData[(d*result_size)+i]));
								__m256d eval_1 = _mm256_load_pd(&(ptrData[(d*result_size)+i+4]));
								__m256d eval_2 = _mm256_load_pd(&(ptrData[(d*result_size)+i+8]));
								__m256d eval_3 = _mm256_load_pd(&(ptrData[(d*result_size)+i+12]));
								__m256d eval_4 = _mm256_load_pd(&(ptrData[(d*result_size)+i+16]));
								__m256d eval_5 = _mm256_load_pd(&(ptrData[(d*result_size)+i+20]));

								__m256d level = _mm256_broadcast_sd(&(ptrLevel[(j*dims)+d]));
								__m256d index = _mm256_broadcast_sd(&(ptrIndex[(j*dims)+d]));

								eval_0 = _mm256_mul_pd(eval_0, level);
								eval_1 = _mm256_mul_pd(eval_1, level);
								eval_2 = _mm256_mul_pd(eval_2, level);
								eval_3 = _mm256_mul_pd(eval_3, level);
								eval_4 = _mm256_mul_pd(eval_4, level);
								eval_5 = _mm256_mul_pd(eval_5, level);

								eval_0 = _mm256_sub_pd(eval_0, index);
								eval_1 = _mm256_sub_pd(eval_1, index);
								eval_2 = _mm256_sub_pd(eval_2, index);
								eval_3 = _mm256_sub_pd(eval_3, index);
								eval_4 = _mm256_sub_pd(eval_4, index);
								eval_5 = _mm256_sub_pd(eval_5, index);

								eval_0 = _mm256_abs_pd(eval_0);
								eval_1 = _mm256_abs_pd(eval_1);
								eval_2 = _mm256_abs_pd(eval_2);
								eval_3 = _mm256_abs_pd(eval_3);
								eval_4 = _mm256_abs_pd(eval_4);
								eval_5 = _mm256_abs_pd(eval_5);

								eval_0 = _mm256_sub_pd(one, eval_0);
								eval_1 = _mm256_sub_pd(one, eval_1);
								eval_2 = _mm256_sub_pd(one, eval_2);
								eval_3 = _mm256_sub_pd(one, eval_3);
								eval_4 = _mm256_sub_pd(one, eval_4);
								eval_5 = _mm256_sub_pd(one, eval_5);

								eval_0 = _mm256_max_pd(zero, eval_0);
								eval_1 = _mm256_max_pd(zero, eval_1);
								eval_2 = _mm256_max_pd(zero, eval_2);
								eval_3 = _mm256_max_pd(zero, eval_3);
								eval_4 = _mm256_max_pd(zero, eval_4);
								eval_5 = _mm256_max_pd(zero, eval_5);

								support_0 = _mm256_mul_pd(support_0, eval_0);
								support_1 = _mm256_mul_pd(support_1, eval_1);
								support_2 = _mm256_mul_pd(support_2, eval_2);
								support_3 = _mm256_mul_pd(support_3, eval_3);
								support_4 = _mm256_mul_pd(support_4, eval_4);
								support_5 = _mm256_mul_pd(support_5, eval_5);
							}

							__m256d res_0 = _mm256_load_pd(&(ptrResult[i]));
							__m256d res_1 = _mm256_load_pd(&(ptrResult[i+4]));
							__m256d res_2 = _mm256_load_pd(&(ptrResult[i+8]));
							__m256d res_3 = _mm256_load_pd(&(ptrResult[i+12]));
							__m256d res_4 = _mm256_load_pd(&(ptrResult[i+16]));
							__m256d res_5 = _mm256_load_pd(&(ptrResult[i+20]));

							res_0 = _mm256_add_pd(res_0, support_0);
							res_1 = _mm256_add_pd(res_1, support_1);
							res_2 = _mm256_add_pd(res_2, support_2);
							res_3 = _mm256_add_pd(res_3, support_3);
							res_4 = _mm256_add_pd(res_4, support_4);
							res_5 = _mm256_add_pd(res_5, support_5);

							_mm256_store_pd(&(ptrResult[i]), res_0);
							_mm256_store_pd(&(ptrResult[i+4]), res_1);
							_mm256_store_pd(&(ptrResult[i+8]), res_2);
							_mm256_store_pd(&(ptrResult[i+12]), res_3);
							_mm256_store_pd(&(ptrResult[i+16]), res_4);
							_mm256_store_pd(&(ptrResult[i+20]), res_5);
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

							#pragma ivdep
							#pragma vector aligned
							for (size_t d = 0; d < dims; d++)
							{
								double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*result_size)+i]));
								double index_calc = eval - (ptrIndex[(j*dims)+d]);
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
							double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*result_size)+i]));
							double index_calc = eval - (ptrIndex[(j*dims)+d]);
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
#ifdef USEOMP
	}
#endif
}

}
