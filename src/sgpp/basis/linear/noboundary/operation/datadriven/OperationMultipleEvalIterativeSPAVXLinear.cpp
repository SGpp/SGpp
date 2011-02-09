/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPAVXLinear.hpp"
#include "exception/operation_exception.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

#ifdef USEICCINTRINSICS
#ifdef USEAVX
#include <immintrin.h>
#else
#include "common/avxintrin_emu.h"
#endif

union floatAbsMaskAVX
{
   const float f;
   const int i;

   floatAbsMaskAVX() : i(0x7FFFFFFF) {}
};

__declspec(align(32)) const floatAbsMaskAVX absMaskSPAVX;

static const __m256 abs2MaskSPAVX = _mm256_broadcast_ss( &(absMaskSPAVX.f) );

const __m256 _mm256_abs_ps( const __m256& x)
{
       return _mm256_and_ps( abs2MaskSPAVX, x);
}

static const __m256i ldStMaskSPAVX = _mm256_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF);
#endif

#define CHUNKDATAPOINTS_AVX 48 // must be divide-able by 48
#define CHUNKGRIDPOINTS_AVX 12

namespace sg
{

OperationMultipleEvalIterativeSPAVXLinear::OperationMultipleEvalIterativeSPAVXLinear(GridStorage* storage, DataMatrixSP* dataset) : OperationMultipleEvalVectorizedSP(dataset)
{
	this->storage = storage;

	this->level_ = new DataMatrixSP(storage->size(), storage->dim());
	this->index_ = new DataMatrixSP(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myTimer = new SGppStopwatch();
}

OperationMultipleEvalIterativeSPAVXLinear::~OperationMultipleEvalIterativeSPAVXLinear()
{
	delete myTimer;
}

void OperationMultipleEvalIterativeSPAVXLinear::rebuildLevelAndIndex()
{
	delete this->level_;
	delete this->index_;

	this->level_ = new DataMatrixSP(storage->size(), storage->dim());
	this->index_ = new DataMatrixSP(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
}

double OperationMultipleEvalIterativeSPAVXLinear::multTransposeVectorized(DataVectorSP& source, DataVectorSP& result)
{
	size_t source_size = source.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();
    float* ptrSource = source.getPointer();
    float* ptrData = this->dataset_->getPointer();
    float* ptrLevel = this->level_->getPointer();
    float* ptrIndex = this->index_->getPointer();
	float* ptrResult = result.getPointer();

    if (this->dataset_->getNcols() % 48 != 0 || source_size != this->dataset_->getNcols())
    {
    	throw operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

    myTimer->start();

    result.setAll(0.0f);

#ifdef _OPENMP
    #pragma omp parallel
	{
		size_t chunksize = (storageSize/omp_get_num_threads())+1;
    	size_t start = chunksize*omp_get_thread_num();
    	size_t end = std::min<size_t>(start+chunksize, storageSize);
#else
    	size_t start = 0;
    	size_t end = storageSize;
#endif
		for(size_t k = start; k < end; k+=std::min<size_t>((size_t)CHUNKGRIDPOINTS_AVX, (end-k)))
		{
			size_t grid_inc = std::min<size_t>((size_t)CHUNKGRIDPOINTS_AVX, (end-k));

#ifdef USEICCINTRINSICS
			for (size_t i = 0; i < source_size; i+=CHUNKDATAPOINTS_AVX)
			{
				for (size_t j = k; j < k+grid_inc; j++)
				{
					__m256 support_0 = _mm256_load_ps(&(ptrSource[i]));
					__m256 support_1 = _mm256_load_ps(&(ptrSource[i+8]));
					__m256 support_2 = _mm256_load_ps(&(ptrSource[i+16]));
					__m256 support_3 = _mm256_load_ps(&(ptrSource[i+24]));
					__m256 support_4 = _mm256_load_ps(&(ptrSource[i+32]));
					__m256 support_5 = _mm256_load_ps(&(ptrSource[i+40]));

					__m256 one = _mm256_set1_ps(1.0f);
					__m256 zero = _mm256_set1_ps(0.0f);

					for (size_t d = 0; d < dims; d++)
					{
						__m256 eval_0 = _mm256_load_ps(&(ptrData[(d*source_size)+i]));
						__m256 eval_1 = _mm256_load_ps(&(ptrData[(d*source_size)+i+8]));
						__m256 eval_2 = _mm256_load_ps(&(ptrData[(d*source_size)+i+16]));
						__m256 eval_3 = _mm256_load_ps(&(ptrData[(d*source_size)+i+24]));
						__m256 eval_4 = _mm256_load_ps(&(ptrData[(d*source_size)+i+32]));
						__m256 eval_5 = _mm256_load_ps(&(ptrData[(d*source_size)+i+40]));

						__m256 level = _mm256_broadcast_ss(&(ptrLevel[(j*dims)+d]));
						__m256 index = _mm256_broadcast_ss(&(ptrIndex[(j*dims)+d]));

						eval_0 = _mm256_mul_ps(eval_0, level);
						eval_1 = _mm256_mul_ps(eval_1, level);
						eval_2 = _mm256_mul_ps(eval_2, level);
						eval_3 = _mm256_mul_ps(eval_3, level);
						eval_4 = _mm256_mul_ps(eval_4, level);
						eval_5 = _mm256_mul_ps(eval_5, level);

						eval_0 = _mm256_sub_ps(eval_0, index);
						eval_1 = _mm256_sub_ps(eval_1, index);
						eval_2 = _mm256_sub_ps(eval_2, index);
						eval_3 = _mm256_sub_ps(eval_3, index);
						eval_4 = _mm256_sub_ps(eval_4, index);
						eval_5 = _mm256_sub_ps(eval_5, index);

						eval_0 = _mm256_abs_ps(eval_0);
						eval_1 = _mm256_abs_ps(eval_1);
						eval_2 = _mm256_abs_ps(eval_2);
						eval_3 = _mm256_abs_ps(eval_3);
						eval_4 = _mm256_abs_ps(eval_4);
						eval_5 = _mm256_abs_ps(eval_5);

						eval_0 = _mm256_sub_ps(one, eval_0);
						eval_1 = _mm256_sub_ps(one, eval_1);
						eval_2 = _mm256_sub_ps(one, eval_2);
						eval_3 = _mm256_sub_ps(one, eval_3);
						eval_4 = _mm256_sub_ps(one, eval_4);
						eval_5 = _mm256_sub_ps(one, eval_5);

						eval_0 = _mm256_max_ps(zero, eval_0);
						eval_1 = _mm256_max_ps(zero, eval_1);
						eval_2 = _mm256_max_ps(zero, eval_2);
						eval_3 = _mm256_max_ps(zero, eval_3);
						eval_4 = _mm256_max_ps(zero, eval_4);
						eval_5 = _mm256_max_ps(zero, eval_5);

						support_0 = _mm256_mul_ps(support_0, eval_0);
						support_1 = _mm256_mul_ps(support_1, eval_1);
						support_2 = _mm256_mul_ps(support_2, eval_2);
						support_3 = _mm256_mul_ps(support_3, eval_3);
						support_4 = _mm256_mul_ps(support_4, eval_4);
						support_5 = _mm256_mul_ps(support_5, eval_5);
					}

					__m256 res_0 = _mm256_maskload_ps(&(ptrResult[j]), ldStMaskSPAVX);

					support_0 = _mm256_add_ps(support_0, support_1);
					support_2 = _mm256_add_ps(support_2, support_3);
					support_4 = _mm256_add_ps(support_4, support_5);
					support_0 = _mm256_add_ps(support_0, support_2);
					support_0 = _mm256_add_ps(support_0, support_4);

					support_0 = _mm256_hadd_ps(support_0, support_0);
					__m256 tmp = _mm256_permute2f128_ps(support_0, support_0, 0x81);
					support_0 = _mm256_add_ps(support_0, tmp);
					support_0 = _mm256_hadd_ps(support_0, support_0);
					res_0 = _mm256_add_ps(res_0, support_0);

					_mm256_maskstore_ps(&(ptrResult[j]), ldStMaskSPAVX, res_0);
				}
			}
#else
			for (size_t i = 0; i < source_size; i++)
			{
				for (size_t j = k; j < k+grid_inc; j++)
				{
					float curSupport = ptrSource[i];

#ifdef __ICC
					#pragma ivdep
					#pragma vector aligned
#endif
					for (size_t d = 0; d < dims; d++)
					{
						float eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*source_size)+i]));
						float index_calc = eval - (ptrIndex[(j*dims)+d]);
						float abs = fabs(index_calc);
						float last = 1.0f - abs;
						float localSupport = std::max<float>(last, 0.0f);
						curSupport *= localSupport;
					}

					ptrResult[j] += curSupport;
				}
			}
#endif
		}
#ifdef _OPENMP
	}
#endif

	return myTimer->stop();
}

double OperationMultipleEvalIterativeSPAVXLinear::multVectorized(DataVectorSP& alpha, DataVectorSP& result)
{
	size_t result_size = result.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();
    float* ptrAlpha = alpha.getPointer();
    float* ptrData = this->dataset_->getPointer();
    float* ptrResult = result.getPointer();
    float* ptrLevel = this->level_->getPointer();
    float* ptrIndex = this->index_->getPointer();

    if (this->dataset_->getNcols() % 48 != 0 || result_size != this->dataset_->getNcols())
    {
    	throw operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

    myTimer->start();

#ifdef _OPENMP
    #pragma omp parallel
	{
		size_t chunksize = (result_size/omp_get_num_threads())+1;
		// assure that every subarray is 32-byte aligned
		if (chunksize % 48 != 0)
		{
			size_t remainder = chunksize % 48;
			size_t patch = 48 - remainder;
			chunksize += patch;
		}
    	size_t start = chunksize*omp_get_thread_num();
    	size_t end = std::min<size_t>(start+chunksize, result_size);
#else
    	size_t start = 0;
    	size_t end = result_size;
#endif
		for(size_t c = start; c < end; c+=std::min<size_t>((size_t)CHUNKDATAPOINTS_AVX, (end-c)))
		{
			size_t data_end = std::min<size_t>((size_t)CHUNKDATAPOINTS_AVX+c, end);

#ifdef __ICC
			#pragma ivdep
			#pragma vector aligned
#endif
			for (size_t i = c; i < data_end; i++)
			{
				ptrResult[i] = 0.0f;
			}

			for (size_t m = 0; m < storageSize; m+=std::min<size_t>((size_t)CHUNKGRIDPOINTS_AVX, (storageSize-m)))
			{
#ifdef USEICCINTRINSICS
				size_t grid_inc = std::min<size_t>((size_t)CHUNKGRIDPOINTS_AVX, (storageSize-m));

				for (size_t i = c; i < c+CHUNKDATAPOINTS_AVX; i+=48)
				{
					for (size_t j = m; j < m+grid_inc; j++)
					{
						__m256 support_0 = _mm256_broadcast_ss(&(ptrAlpha[j]));
						__m256 support_1 = _mm256_broadcast_ss(&(ptrAlpha[j]));
						__m256 support_2 = _mm256_broadcast_ss(&(ptrAlpha[j]));
						__m256 support_3 = _mm256_broadcast_ss(&(ptrAlpha[j]));
						__m256 support_4 = _mm256_broadcast_ss(&(ptrAlpha[j]));
						__m256 support_5 = _mm256_broadcast_ss(&(ptrAlpha[j]));

						__m256 one = _mm256_set1_ps(1.0f);
						__m256 zero = _mm256_set1_ps(0.0f);

						for (size_t d = 0; d < dims; d++)
						{
							__m256 eval_0 = _mm256_load_ps(&(ptrData[(d*result_size)+i]));
							__m256 eval_1 = _mm256_load_ps(&(ptrData[(d*result_size)+i+8]));
							__m256 eval_2 = _mm256_load_ps(&(ptrData[(d*result_size)+i+16]));
							__m256 eval_3 = _mm256_load_ps(&(ptrData[(d*result_size)+i+24]));
							__m256 eval_4 = _mm256_load_ps(&(ptrData[(d*result_size)+i+32]));
							__m256 eval_5 = _mm256_load_ps(&(ptrData[(d*result_size)+i+40]));

							__m256 level = _mm256_broadcast_ss(&(ptrLevel[(j*dims)+d]));
							__m256 index = _mm256_broadcast_ss(&(ptrIndex[(j*dims)+d]));

							eval_0 = _mm256_mul_ps(eval_0, level);
							eval_1 = _mm256_mul_ps(eval_1, level);
							eval_2 = _mm256_mul_ps(eval_2, level);
							eval_3 = _mm256_mul_ps(eval_3, level);
							eval_4 = _mm256_mul_ps(eval_4, level);
							eval_5 = _mm256_mul_ps(eval_5, level);

							eval_0 = _mm256_sub_ps(eval_0, index);
							eval_1 = _mm256_sub_ps(eval_1, index);
							eval_2 = _mm256_sub_ps(eval_2, index);
							eval_3 = _mm256_sub_ps(eval_3, index);
							eval_4 = _mm256_sub_ps(eval_4, index);
							eval_5 = _mm256_sub_ps(eval_5, index);

							eval_0 = _mm256_abs_ps(eval_0);
							eval_1 = _mm256_abs_ps(eval_1);
							eval_2 = _mm256_abs_ps(eval_2);
							eval_3 = _mm256_abs_ps(eval_3);
							eval_4 = _mm256_abs_ps(eval_4);
							eval_5 = _mm256_abs_ps(eval_5);

							eval_0 = _mm256_sub_ps(one, eval_0);
							eval_1 = _mm256_sub_ps(one, eval_1);
							eval_2 = _mm256_sub_ps(one, eval_2);
							eval_3 = _mm256_sub_ps(one, eval_3);
							eval_4 = _mm256_sub_ps(one, eval_4);
							eval_5 = _mm256_sub_ps(one, eval_5);

							eval_0 = _mm256_max_ps(zero, eval_0);
							eval_1 = _mm256_max_ps(zero, eval_1);
							eval_2 = _mm256_max_ps(zero, eval_2);
							eval_3 = _mm256_max_ps(zero, eval_3);
							eval_4 = _mm256_max_ps(zero, eval_4);
							eval_5 = _mm256_max_ps(zero, eval_5);

							support_0 = _mm256_mul_ps(support_0, eval_0);
							support_1 = _mm256_mul_ps(support_1, eval_1);
							support_2 = _mm256_mul_ps(support_2, eval_2);
							support_3 = _mm256_mul_ps(support_3, eval_3);
							support_4 = _mm256_mul_ps(support_4, eval_4);
							support_5 = _mm256_mul_ps(support_5, eval_5);
						}

						__m256 res_0 = _mm256_load_ps(&(ptrResult[i]));
						__m256 res_1 = _mm256_load_ps(&(ptrResult[i+8]));
						__m256 res_2 = _mm256_load_ps(&(ptrResult[i+16]));
						__m256 res_3 = _mm256_load_ps(&(ptrResult[i+24]));
						__m256 res_4 = _mm256_load_ps(&(ptrResult[i+32]));
						__m256 res_5 = _mm256_load_ps(&(ptrResult[i+40]));

						res_0 = _mm256_add_ps(res_0, support_0);
						res_1 = _mm256_add_ps(res_1, support_1);
						res_2 = _mm256_add_ps(res_2, support_2);
						res_3 = _mm256_add_ps(res_3, support_3);
						res_4 = _mm256_add_ps(res_4, support_4);
						res_5 = _mm256_add_ps(res_5, support_5);

						_mm256_store_ps(&(ptrResult[i]), res_0);
						_mm256_store_ps(&(ptrResult[i+8]), res_1);
						_mm256_store_ps(&(ptrResult[i+16]), res_2);
						_mm256_store_ps(&(ptrResult[i+24]), res_3);
						_mm256_store_ps(&(ptrResult[i+32]), res_4);
						_mm256_store_ps(&(ptrResult[i+40]), res_5);
					}
				}
#else
				size_t grid_end = std::min<size_t>((size_t)CHUNKGRIDPOINTS_AVX+m, storageSize);

				for (size_t i = c; i < data_end; i++)
				{
					for (size_t j = m; j < grid_end; j++)
					{
						float curSupport = ptrAlpha[j];

#ifdef __ICC
						#pragma ivdep
						#pragma vector aligned
#endif
						for (size_t d = 0; d < dims; d++)
						{
							float eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*result_size)+i]));
							float index_calc = eval - (ptrIndex[(j*dims)+d]);
							float abs = fabs(index_calc);
							float last = 1.0f - abs;
							float localSupport = std::max<float>(last, 0.0f);
							curSupport *= localSupport;
						}

						ptrResult[i] += curSupport;
					}
				}
#endif
	        }
		}
#ifdef _OPENMP
	}
#endif

	return myTimer->stop();
}

}
