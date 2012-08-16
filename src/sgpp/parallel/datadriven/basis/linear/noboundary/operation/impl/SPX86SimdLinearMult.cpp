/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/datadriven/basis/linear/noboundary/operation/impl/SPX86SimdLinearMult.h"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/ChunkSizes.h"
#include "base/exception/operation_exception.hpp"

#if defined(__SSE3__) || defined(__AVX__)
#ifdef _WIN32
#include <immintrin.h>
#else
#include <x86intrin.h>
#endif
#endif

#ifdef __USEAVX128__
#undef __AVX__
#endif

namespace sg {
namespace parallel {

SPX86SimdLinearMult::SPX86SimdLinearMult(base::DataMatrixSP *level, base::DataMatrixSP *index, base::DataMatrixSP *dataset, base::DataVectorSP &alpha, base::DataVectorSP &result):
	chunkDataPoints(CHUNKDATAPOINTS_SP_X86),
	chunkGridPoints(CHUNKGRIDPOINTS_SP_X86),
	_level(level),
	_index(index),
	_dataset(dataset),
	_alpha(alpha),
	_result(result)
{
	if (this->_dataset->getNcols() % chunkDataPoints != 0 || _result.getSize() != this->_dataset->getNcols())
	{
		throw sg::base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
	}
}

void SPX86SimdLinearMult::operator ()(size_t start_index_data, size_t end_index_data)
{
	float* ptrLevel = _level->getPointer();
	float* ptrIndex = _index->getPointer();
	float* ptrAlpha = _alpha.getPointer();
	float* ptrData = _dataset->getPointer();
	float* ptrResult = _result.getPointer();
	size_t result_size = _result.getSize();
	size_t dims = _dataset->getNrows();
	size_t storageSize = _alpha.getSize();

	size_t end = end_index_data;

	for(size_t c = start_index_data; c < end; c+=std::min<size_t>(CHUNKDATAPOINTS_SP_X86, (end-c)))
	{
		size_t data_end = std::min<size_t>((size_t)CHUNKDATAPOINTS_SP_X86+c, end);

#ifdef __ICC
#pragma ivdep
#pragma vector aligned
#endif
		for (size_t i = c; i < data_end; i++)
		{
			ptrResult[i] = 0.0f;
		}

		for (size_t m = 0; m < storageSize; m+=std::min<size_t>((size_t)CHUNKGRIDPOINTS_SP_X86, (storageSize-m)))
		{
#if defined(__SSE3__) && !defined(__AVX__)
			size_t grid_inc = std::min<size_t>((size_t)CHUNKGRIDPOINTS_SP_X86, (storageSize-m));

			int imask = 0x7FFFFFFF;
			float* fmask = (float*)&imask;

			for (size_t i = c; i < c+CHUNKDATAPOINTS_SP_X86; i+=24)
			{
				for (size_t j = m; j < m+grid_inc; j++)
				{
					__m128 support_0 = _mm_load1_ps(&(ptrAlpha[j]));
					__m128 support_1 = _mm_load1_ps(&(ptrAlpha[j]));
					__m128 support_2 = _mm_load1_ps(&(ptrAlpha[j]));
					__m128 support_3 = _mm_load1_ps(&(ptrAlpha[j]));
					__m128 support_4 = _mm_load1_ps(&(ptrAlpha[j]));
					__m128 support_5 = _mm_load1_ps(&(ptrAlpha[j]));

					__m128 mask = _mm_set1_ps(*fmask);
					__m128 one = _mm_set1_ps(1.0f);
					__m128 zero = _mm_set1_ps(0.0f);

					for (size_t d = 0; d < dims; d++)
					{
						__m128 eval_0 = _mm_load_ps(&(ptrData[(d*result_size)+i]));
						__m128 eval_1 = _mm_load_ps(&(ptrData[(d*result_size)+i+4]));
						__m128 eval_2 = _mm_load_ps(&(ptrData[(d*result_size)+i+8]));
						__m128 eval_3 = _mm_load_ps(&(ptrData[(d*result_size)+i+12]));
						__m128 eval_4 = _mm_load_ps(&(ptrData[(d*result_size)+i+16]));
						__m128 eval_5 = _mm_load_ps(&(ptrData[(d*result_size)+i+20]));

						__m128 level = _mm_load1_ps(&(ptrLevel[(j*dims)+d]));
						__m128 index = _mm_load1_ps(&(ptrIndex[(j*dims)+d]));
#ifdef __FMA4__
						eval_0 = _mm_msub_ps(eval_0, level, index);
						eval_1 = _mm_msub_ps(eval_1, level, index);
						eval_2 = _mm_msub_ps(eval_2, level, index);
						eval_3 = _mm_msub_ps(eval_3, level, index);
						eval_4 = _mm_msub_ps(eval_4, level, index);
						eval_5 = _mm_msub_ps(eval_5, level, index);
#else
						eval_0 = _mm_sub_ps(_mm_mul_ps(eval_0, level), index);
						eval_1 = _mm_sub_ps(_mm_mul_ps(eval_1, level), index);
						eval_2 = _mm_sub_ps(_mm_mul_ps(eval_2, level), index);
						eval_3 = _mm_sub_ps(_mm_mul_ps(eval_3, level), index);
						eval_4 = _mm_sub_ps(_mm_mul_ps(eval_4, level), index);
						eval_5 = _mm_sub_ps(_mm_mul_ps(eval_5, level), index);
#endif
						eval_0 = _mm_and_ps(mask, eval_0);
						eval_1 = _mm_and_ps(mask, eval_1);
						eval_2 = _mm_and_ps(mask, eval_2);
						eval_3 = _mm_and_ps(mask, eval_3);
						eval_4 = _mm_and_ps(mask, eval_4);
						eval_5 = _mm_and_ps(mask, eval_5);

						eval_0 = _mm_sub_ps(one, eval_0);
						eval_1 = _mm_sub_ps(one, eval_1);
						eval_2 = _mm_sub_ps(one, eval_2);
						eval_3 = _mm_sub_ps(one, eval_3);
						eval_4 = _mm_sub_ps(one, eval_4);
						eval_5 = _mm_sub_ps(one, eval_5);

						eval_0 = _mm_max_ps(zero, eval_0);
						eval_1 = _mm_max_ps(zero, eval_1);
						eval_2 = _mm_max_ps(zero, eval_2);
						eval_3 = _mm_max_ps(zero, eval_3);
						eval_4 = _mm_max_ps(zero, eval_4);
						eval_5 = _mm_max_ps(zero, eval_5);

						support_0 = _mm_mul_ps(support_0, eval_0);
						support_1 = _mm_mul_ps(support_1, eval_1);
						support_2 = _mm_mul_ps(support_2, eval_2);
						support_3 = _mm_mul_ps(support_3, eval_3);
						support_4 = _mm_mul_ps(support_4, eval_4);
						support_5 = _mm_mul_ps(support_5, eval_5);
					}

					__m128 res_0 = _mm_load_ps(&(ptrResult[i]));
					__m128 res_1 = _mm_load_ps(&(ptrResult[i+4]));
					__m128 res_2 = _mm_load_ps(&(ptrResult[i+8]));
					__m128 res_3 = _mm_load_ps(&(ptrResult[i+12]));
					__m128 res_4 = _mm_load_ps(&(ptrResult[i+16]));
					__m128 res_5 = _mm_load_ps(&(ptrResult[i+20]));

					res_0 = _mm_add_ps(res_0, support_0);
					res_1 = _mm_add_ps(res_1, support_1);
					res_2 = _mm_add_ps(res_2, support_2);
					res_3 = _mm_add_ps(res_3, support_3);
					res_4 = _mm_add_ps(res_4, support_4);
					res_5 = _mm_add_ps(res_5, support_5);

					_mm_store_ps(&(ptrResult[i]), res_0);
					_mm_store_ps(&(ptrResult[i+4]), res_1);
					_mm_store_ps(&(ptrResult[i+8]), res_2);
					_mm_store_ps(&(ptrResult[i+12]), res_3);
					_mm_store_ps(&(ptrResult[i+16]), res_4);
					_mm_store_ps(&(ptrResult[i+20]), res_5);
				}
			}
#endif
#if defined(__SSE3__) && defined(__AVX__)
			size_t grid_inc = std::min<size_t>((size_t)CHUNKGRIDPOINTS_SP_X86, (storageSize-m));

			int imask = 0x7FFFFFFF;
			float* fmask = (float*)&imask;

			for (size_t i = c; i < c+CHUNKDATAPOINTS_SP_X86; i+=48)
			{
				for (size_t j = m; j < m+grid_inc; j++)
				{
					__m256 support_0 = _mm256_broadcast_ss(&(ptrAlpha[j]));
					__m256 support_1 = _mm256_broadcast_ss(&(ptrAlpha[j]));
					__m256 support_2 = _mm256_broadcast_ss(&(ptrAlpha[j]));
					__m256 support_3 = _mm256_broadcast_ss(&(ptrAlpha[j]));
					__m256 support_4 = _mm256_broadcast_ss(&(ptrAlpha[j]));
					__m256 support_5 = _mm256_broadcast_ss(&(ptrAlpha[j]));

					__m256 mask = _mm256_set1_ps(*fmask);
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
#ifdef __FMA4__
						eval_0 = _mm256_msub_ps(eval_0, level, index);
						eval_1 = _mm256_msub_ps(eval_1, level, index);
						eval_2 = _mm256_msub_ps(eval_2, level, index);
						eval_3 = _mm256_msub_ps(eval_3, level, index);
						eval_4 = _mm256_msub_ps(eval_4, level, index);
						eval_5 = _mm256_msub_ps(eval_5, level, index);
#else
						eval_0 = _mm256_sub_ps(_mm256_mul_ps(eval_0, level), index);
						eval_1 = _mm256_sub_ps(_mm256_mul_ps(eval_1, level), index);
						eval_2 = _mm256_sub_ps(_mm256_mul_ps(eval_2, level), index);
						eval_3 = _mm256_sub_ps(_mm256_mul_ps(eval_3, level), index);
						eval_4 = _mm256_sub_ps(_mm256_mul_ps(eval_4, level), index);
						eval_5 = _mm256_sub_ps(_mm256_mul_ps(eval_5, level), index);
#endif
						eval_0 = _mm256_and_ps(mask, eval_0);
						eval_1 = _mm256_and_ps(mask, eval_1);
						eval_2 = _mm256_and_ps(mask, eval_2);
						eval_3 = _mm256_and_ps(mask, eval_3);
						eval_4 = _mm256_and_ps(mask, eval_4);
						eval_5 = _mm256_and_ps(mask, eval_5);

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
#endif
#if !defined(__SSE3__) && !defined(__AVX__)
			size_t grid_end = std::min<size_t>((size_t)CHUNKGRIDPOINTS_SP_X86+m, storageSize);

			for (size_t i = c; i < data_end; i++)
			{
				for (size_t j = m; j < grid_end; j++)
				{
					float curSupport = ptrAlpha[j];

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
}


}
}
