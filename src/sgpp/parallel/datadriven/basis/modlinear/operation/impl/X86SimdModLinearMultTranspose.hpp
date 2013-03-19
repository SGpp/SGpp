/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef X86SIMDMODLINEARMULTTRANSPOSE_HPP
#define X86SIMDMODLINEARMULTTRANSPOSE_HPP

#include "base/grid/GridStorage.hpp"

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

class X86SimdModLinearMultTranspose
{
public:
	static inline size_t getChunkGridPoints(){return 12;}
	static inline size_t getChunkDataPoints(){return 24;}
	static inline void multTranspose(
			sg::base::DataMatrix* level,
			sg::base::DataMatrix* index,
			sg::base::DataMatrix* dataset,
			sg::base::DataVector& source,
			sg::base::DataVector& result,
			size_t start_index_grid,
			size_t end_index_grid,
			size_t start_index_data,
			size_t end_index_data){
		double* ptrLevel = level->getPointer();
		double* ptrIndex = index->getPointer();
		double* ptrSource = source.getPointer();
		double* ptrData = dataset->getPointer();
		double* ptrResult = result.getPointer();

		size_t source_size = source.getSize();
		size_t dims = dataset->getNrows();

		size_t end = end_index_grid;

		for(size_t k = start_index_grid; k < end; k+=std::min<size_t>((size_t)getChunkGridPoints(), (end-k)))
		{
			size_t grid_inc = std::min<size_t>((size_t)getChunkGridPoints(), (end-k));
	#if defined(__SSE3__) && !defined(__AVX__)
			long long imask = 0x7FFFFFFFFFFFFFFF;
			double* fmask = (double*)&imask;

			for (size_t i = start_index_data; i < end_index_data; i+=12)
			{
				for (size_t j = k; j < k+grid_inc; j++)
				{
					__m128d support_0 = _mm_load_pd(&(ptrSource[i]));
					__m128d support_1 = _mm_load_pd(&(ptrSource[i+2]));
					__m128d support_2 = _mm_load_pd(&(ptrSource[i+4]));
					__m128d support_3 = _mm_load_pd(&(ptrSource[i+6]));
					__m128d support_4 = _mm_load_pd(&(ptrSource[i+8]));
					__m128d support_5 = _mm_load_pd(&(ptrSource[i+10]));

					__m128d one = _mm_set1_pd(1.0);
					__m128d two = _mm_set1_pd(2.0);
					__m128d zero = _mm_set1_pd(0.0);

					for (size_t d = 0; d < dims; d++)
					{
						// special case for level 1
						if (ptrLevel[(j*dims)+d] == 2.0)
						{
							// Nothing (multiply by one)
						}
						// most left basis function on every level
						else if (ptrIndex[(j*dims)+d] == 1.0)
						{
							__m128d eval_0 = _mm_load_pd(&(ptrData[(d*source_size)+i]));
							__m128d eval_1 = _mm_load_pd(&(ptrData[(d*source_size)+i+2]));
							__m128d eval_2 = _mm_load_pd(&(ptrData[(d*source_size)+i+4]));
							__m128d eval_3 = _mm_load_pd(&(ptrData[(d*source_size)+i+6]));
							__m128d eval_4 = _mm_load_pd(&(ptrData[(d*source_size)+i+8]));
							__m128d eval_5 = _mm_load_pd(&(ptrData[(d*source_size)+i+10]));

							__m128d level = _mm_loaddup_pd(&(ptrLevel[(j*dims)+d]));

							eval_0 = _mm_mul_pd(eval_0, level);
							eval_1 = _mm_mul_pd(eval_1, level);
							eval_2 = _mm_mul_pd(eval_2, level);
							eval_3 = _mm_mul_pd(eval_3, level);
							eval_4 = _mm_mul_pd(eval_4, level);
							eval_5 = _mm_mul_pd(eval_5, level);

							eval_0 = _mm_sub_pd(two, eval_0);
							eval_1 = _mm_sub_pd(two, eval_1);
							eval_2 = _mm_sub_pd(two, eval_2);
							eval_3 = _mm_sub_pd(two, eval_3);
							eval_4 = _mm_sub_pd(two, eval_4);
							eval_5 = _mm_sub_pd(two, eval_5);

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
						// most right basis function on every level
						else if (ptrIndex[(j*dims)+d] == (ptrLevel[(j*dims)+d] - 1.0))
						{
							__m128d eval_0 = _mm_load_pd(&(ptrData[(d*source_size)+i]));
							__m128d eval_1 = _mm_load_pd(&(ptrData[(d*source_size)+i+2]));
							__m128d eval_2 = _mm_load_pd(&(ptrData[(d*source_size)+i+4]));
							__m128d eval_3 = _mm_load_pd(&(ptrData[(d*source_size)+i+6]));
							__m128d eval_4 = _mm_load_pd(&(ptrData[(d*source_size)+i+8]));
							__m128d eval_5 = _mm_load_pd(&(ptrData[(d*source_size)+i+10]));

							__m128d level = _mm_loaddup_pd(&(ptrLevel[(j*dims)+d]));
							__m128d index = _mm_loaddup_pd(&(ptrIndex[(j*dims)+d]));
	#ifdef __FMA4__
							eval_0 = _mm_msub_pd(eval_0, level, index);
							eval_1 = _mm_msub_pd(eval_1, level, index);
							eval_2 = _mm_msub_pd(eval_2, level, index);
							eval_3 = _mm_msub_pd(eval_3, level, index);
							eval_4 = _mm_msub_pd(eval_4, level, index);
							eval_5 = _mm_msub_pd(eval_5, level, index);
	#else
							eval_0 = _mm_sub_pd(_mm_mul_pd(eval_0, level), index);
							eval_1 = _mm_sub_pd(_mm_mul_pd(eval_1, level), index);
							eval_2 = _mm_sub_pd(_mm_mul_pd(eval_2, level), index);
							eval_3 = _mm_sub_pd(_mm_mul_pd(eval_3, level), index);
							eval_4 = _mm_sub_pd(_mm_mul_pd(eval_4, level), index);
							eval_5 = _mm_sub_pd(_mm_mul_pd(eval_5, level), index);
	#endif
							eval_0 = _mm_add_pd(one, eval_0);
							eval_1 = _mm_add_pd(one, eval_1);
							eval_2 = _mm_add_pd(one, eval_2);
							eval_3 = _mm_add_pd(one, eval_3);
							eval_4 = _mm_add_pd(one, eval_4);
							eval_5 = _mm_add_pd(one, eval_5);

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
						// all other basis functions
						else
						{
							__m128d eval_0 = _mm_load_pd(&(ptrData[(d*source_size)+i]));
							__m128d eval_1 = _mm_load_pd(&(ptrData[(d*source_size)+i+2]));
							__m128d eval_2 = _mm_load_pd(&(ptrData[(d*source_size)+i+4]));
							__m128d eval_3 = _mm_load_pd(&(ptrData[(d*source_size)+i+6]));
							__m128d eval_4 = _mm_load_pd(&(ptrData[(d*source_size)+i+8]));
							__m128d eval_5 = _mm_load_pd(&(ptrData[(d*source_size)+i+10]));

							__m128d level = _mm_loaddup_pd(&(ptrLevel[(j*dims)+d]));
							__m128d index = _mm_loaddup_pd(&(ptrIndex[(j*dims)+d]));
	#ifdef __FMA4__
							eval_0 = _mm_msub_pd(eval_0, level, index);
							eval_1 = _mm_msub_pd(eval_1, level, index);
							eval_2 = _mm_msub_pd(eval_2, level, index);
							eval_3 = _mm_msub_pd(eval_3, level, index);
							eval_4 = _mm_msub_pd(eval_4, level, index);
							eval_5 = _mm_msub_pd(eval_5, level, index);
	#else
							eval_0 = _mm_sub_pd(_mm_mul_pd(eval_0, level), index);
							eval_1 = _mm_sub_pd(_mm_mul_pd(eval_1, level), index);
							eval_2 = _mm_sub_pd(_mm_mul_pd(eval_2, level), index);
							eval_3 = _mm_sub_pd(_mm_mul_pd(eval_3, level), index);
							eval_4 = _mm_sub_pd(_mm_mul_pd(eval_4, level), index);
							eval_5 = _mm_sub_pd(_mm_mul_pd(eval_5, level), index);
	#endif
							__m128d mask = _mm_set1_pd(*fmask);

							eval_0 = _mm_and_pd(mask, eval_0);
							eval_1 = _mm_and_pd(mask, eval_1);
							eval_2 = _mm_and_pd(mask, eval_2);
							eval_3 = _mm_and_pd(mask, eval_3);
							eval_4 = _mm_and_pd(mask, eval_4);
							eval_5 = _mm_and_pd(mask, eval_5);

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
	#endif
	#if defined(__SSE3__) && defined(__AVX__)
			long long imask = 0x7FFFFFFFFFFFFFFF;
			double* fmask = (double*)&imask;

			for (size_t i = start_index_data; i < end_index_data; i+=24)
			{
				for (size_t j = k; j < k+grid_inc; j++)
				{
					__m256d support_0 = _mm256_load_pd(&(ptrSource[i]));
					__m256d support_1 = _mm256_load_pd(&(ptrSource[i+4]));
					__m256d support_2 = _mm256_load_pd(&(ptrSource[i+8]));
					__m256d support_3 = _mm256_load_pd(&(ptrSource[i+12]));
					__m256d support_4 = _mm256_load_pd(&(ptrSource[i+16]));
					__m256d support_5 = _mm256_load_pd(&(ptrSource[i+20]));

					__m256d one = _mm256_set1_pd(1.0);
					__m256d two = _mm256_set1_pd(2.0);
					__m256d zero = _mm256_set1_pd(0.0);

					for (size_t d = 0; d < dims; d++)
					{
						// special case for level 1
						if (ptrLevel[(j*dims)+d] == 2.0)
						{
							// Nothing (multiply by one)
						}
						// most left basis function on every level
						else if (ptrIndex[(j*dims)+d] == 1.0)
						{
							__m256d eval_0 = _mm256_load_pd(&(ptrData[(d*source_size)+i]));
							__m256d eval_1 = _mm256_load_pd(&(ptrData[(d*source_size)+i+4]));
							__m256d eval_2 = _mm256_load_pd(&(ptrData[(d*source_size)+i+8]));
							__m256d eval_3 = _mm256_load_pd(&(ptrData[(d*source_size)+i+12]));
							__m256d eval_4 = _mm256_load_pd(&(ptrData[(d*source_size)+i+16]));
							__m256d eval_5 = _mm256_load_pd(&(ptrData[(d*source_size)+i+20]));

							__m256d level = _mm256_broadcast_sd(&(ptrLevel[(j*dims)+d]));

							eval_0 = _mm256_mul_pd(eval_0, level);
							eval_1 = _mm256_mul_pd(eval_1, level);
							eval_2 = _mm256_mul_pd(eval_2, level);
							eval_3 = _mm256_mul_pd(eval_3, level);
							eval_4 = _mm256_mul_pd(eval_4, level);
							eval_5 = _mm256_mul_pd(eval_5, level);

							eval_0 = _mm256_sub_pd(two, eval_0);
							eval_1 = _mm256_sub_pd(two, eval_1);
							eval_2 = _mm256_sub_pd(two, eval_2);
							eval_3 = _mm256_sub_pd(two, eval_3);
							eval_4 = _mm256_sub_pd(two, eval_4);
							eval_5 = _mm256_sub_pd(two, eval_5);

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
						// most right basis function on every level
						else if (ptrIndex[(j*dims)+d] == (ptrLevel[(j*dims)+d] - 1.0f))
						{
							__m256d eval_0 = _mm256_load_pd(&(ptrData[(d*source_size)+i]));
							__m256d eval_1 = _mm256_load_pd(&(ptrData[(d*source_size)+i+4]));
							__m256d eval_2 = _mm256_load_pd(&(ptrData[(d*source_size)+i+8]));
							__m256d eval_3 = _mm256_load_pd(&(ptrData[(d*source_size)+i+12]));
							__m256d eval_4 = _mm256_load_pd(&(ptrData[(d*source_size)+i+16]));
							__m256d eval_5 = _mm256_load_pd(&(ptrData[(d*source_size)+i+20]));

							__m256d level = _mm256_broadcast_sd(&(ptrLevel[(j*dims)+d]));
							__m256d index = _mm256_broadcast_sd(&(ptrIndex[(j*dims)+d]));
	#ifdef __FMA4__
							eval_0 = _mm256_msub_pd(eval_0, level, index);
							eval_1 = _mm256_msub_pd(eval_1, level, index);
							eval_2 = _mm256_msub_pd(eval_2, level, index);
							eval_3 = _mm256_msub_pd(eval_3, level, index);
							eval_4 = _mm256_msub_pd(eval_4, level, index);
							eval_5 = _mm256_msub_pd(eval_5, level, index);
	#else
							eval_0 = _mm256_sub_pd(_mm256_mul_pd(eval_0, level), index);
							eval_1 = _mm256_sub_pd(_mm256_mul_pd(eval_1, level), index);
							eval_2 = _mm256_sub_pd(_mm256_mul_pd(eval_2, level), index);
							eval_3 = _mm256_sub_pd(_mm256_mul_pd(eval_3, level), index);
							eval_4 = _mm256_sub_pd(_mm256_mul_pd(eval_4, level), index);
							eval_5 = _mm256_sub_pd(_mm256_mul_pd(eval_5, level), index);
	#endif
							eval_0 = _mm256_add_pd(one, eval_0);
							eval_1 = _mm256_add_pd(one, eval_1);
							eval_2 = _mm256_add_pd(one, eval_2);
							eval_3 = _mm256_add_pd(one, eval_3);
							eval_4 = _mm256_add_pd(one, eval_4);
							eval_5 = _mm256_add_pd(one, eval_5);

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
						// all other basis functions
						else
						{
							__m256d eval_0 = _mm256_load_pd(&(ptrData[(d*source_size)+i]));
							__m256d eval_1 = _mm256_load_pd(&(ptrData[(d*source_size)+i+4]));
							__m256d eval_2 = _mm256_load_pd(&(ptrData[(d*source_size)+i+8]));
							__m256d eval_3 = _mm256_load_pd(&(ptrData[(d*source_size)+i+12]));
							__m256d eval_4 = _mm256_load_pd(&(ptrData[(d*source_size)+i+16]));
							__m256d eval_5 = _mm256_load_pd(&(ptrData[(d*source_size)+i+20]));

							__m256d level = _mm256_broadcast_sd(&(ptrLevel[(j*dims)+d]));
							__m256d index = _mm256_broadcast_sd(&(ptrIndex[(j*dims)+d]));
	#ifdef __FMA4__
							eval_0 = _mm256_msub_pd(eval_0, level, index);
							eval_1 = _mm256_msub_pd(eval_1, level, index);
							eval_2 = _mm256_msub_pd(eval_2, level, index);
							eval_3 = _mm256_msub_pd(eval_3, level, index);
							eval_4 = _mm256_msub_pd(eval_4, level, index);
							eval_5 = _mm256_msub_pd(eval_5, level, index);
	#else
							eval_0 = _mm256_sub_pd(_mm256_mul_pd(eval_0, level), index);
							eval_1 = _mm256_sub_pd(_mm256_mul_pd(eval_1, level), index);
							eval_2 = _mm256_sub_pd(_mm256_mul_pd(eval_2, level), index);
							eval_3 = _mm256_sub_pd(_mm256_mul_pd(eval_3, level), index);
							eval_4 = _mm256_sub_pd(_mm256_mul_pd(eval_4, level), index);
							eval_5 = _mm256_sub_pd(_mm256_mul_pd(eval_5, level), index);
	#endif
							__m256d mask = _mm256_broadcast_sd(fmask);

							eval_0 = _mm256_and_pd(mask, eval_0);
							eval_1 = _mm256_and_pd(mask, eval_1);
							eval_2 = _mm256_and_pd(mask, eval_2);
							eval_3 = _mm256_and_pd(mask, eval_3);
							eval_4 = _mm256_and_pd(mask, eval_4);
							eval_5 = _mm256_and_pd(mask, eval_5);

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
					}

					const __m256i ldStMaskAVX = _mm256_set_epi64x(0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xFFFFFFFFFFFFFFFF);

					support_0 = _mm256_add_pd(support_0, support_1);
					support_2 = _mm256_add_pd(support_2, support_3);
					support_4 = _mm256_add_pd(support_4, support_5);
					support_0 = _mm256_add_pd(support_0, support_2);
					support_0 = _mm256_add_pd(support_0, support_4);

					support_0 = _mm256_hadd_pd(support_0, support_0);
					__m256d tmp = _mm256_permute2f128_pd(support_0, support_0, 0x81);
					support_0 = _mm256_add_pd(support_0, tmp);

	// Workaround: bug with maskload in GCC (4.6.1)
	#ifdef __ICC
					__m256d res_0 = _mm256_maskload_pd(&(ptrResult[j]), ldStMaskAVX);
					res_0 = _mm256_add_pd(res_0, support_0);
					_mm256_maskstore_pd(&(ptrResult[j]), ldStMaskAVX, res_0);
	#else
					double tmp_reduce;
					_mm256_maskstore_pd(&(tmp_reduce), ldStMaskAVX, support_0);
					ptrResult[j] += tmp_reduce;
	#endif
				}
			}
	#endif
	#if !defined(__SSE3__) && !defined(__AVX__)
			for (size_t i = start_index_data; i < end_index_data; i++)
			{
				for (size_t j = k; j < k+grid_inc; j++)
				{
					double curSupport = ptrSource[i];

					for (size_t d = 0; d < dims; d++)
					{
						if (ptrLevel[(j*dims)+d] == 2.0)
						{
							// nothing to do (mult with 1)
						}
						else if (ptrIndex[(j*dims)+d] == 1.0)
						{
							double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*source_size)+i]));
							eval = 2.0 - eval;
							double localSupport = std::max<double>(eval, 0.0);
							curSupport *= localSupport;
						}
						else if (ptrIndex[(j*dims)+d] == (ptrLevel[(j*dims)+d] - 1.0))
						{
							double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*source_size)+i]));
							double index_calc = eval - (ptrIndex[(j*dims)+d]);
							double last = 1.0 + index_calc;
							double localSupport = std::max<double>(last, 0.0);
							curSupport *= localSupport;
						}
						else
						{
							double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*source_size)+i]));
							double index_calc = eval - (ptrIndex[(j*dims)+d]);
							double abs = fabs(index_calc);
							double last = 1.0 - abs;
							double localSupport = std::max<double>(last, 0.0);
							curSupport *= localSupport;
						}
					}

					ptrResult[j] += curSupport;
				}
			}
	#endif
		}
	}
};

}
}

#endif // X86SIMDMODLINEARMULTTRANSPOSE_HPP
