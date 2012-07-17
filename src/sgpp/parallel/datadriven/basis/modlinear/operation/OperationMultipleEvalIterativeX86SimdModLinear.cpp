/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeX86SimdModLinear.hpp"
#include "base/exception/operation_exception.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

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

#define CHUNKDATAPOINTS_X86 24 // must be divide-able by 24
#define CHUNKGRIDPOINTS_x86 12

namespace sg
{
namespace parallel
{

OperationMultipleEvalIterativeX86SimdModLinear::OperationMultipleEvalIterativeX86SimdModLinear(sg::base::GridStorage* storage, sg::base::DataMatrix* dataset) : sg::parallel::OperationMultipleEvalVectorized(dataset)
{
	this->storage = storage;

	this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
	this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myTimer = new sg::base::SGppStopwatch();
}

OperationMultipleEvalIterativeX86SimdModLinear::~OperationMultipleEvalIterativeX86SimdModLinear()
{
	delete myTimer;
}

void OperationMultipleEvalIterativeX86SimdModLinear::rebuildLevelAndIndex()
{
	delete this->level_;
	delete this->index_;

	this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
	this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
}

double OperationMultipleEvalIterativeX86SimdModLinear::multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result)
{
	size_t source_size = source.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();
    double* ptrSource = source.getPointer();
    double* ptrData = this->dataset_->getPointer();
    double* ptrLevel = this->level_->getPointer();
    double* ptrIndex = this->index_->getPointer();
    double* ptrResult = result.getPointer();

    if (this->dataset_->getNcols() % CHUNKDATAPOINTS_X86 != 0 || source_size != this->dataset_->getNcols())
    {
    	throw sg::base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

    myTimer->start();

    result.setAll(0.0);

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
		for(size_t k = start; k < end; k+=std::min<size_t>((size_t)CHUNKGRIDPOINTS_x86, (end-k)))
		{
			size_t grid_inc = std::min<size_t>((size_t)CHUNKGRIDPOINTS_x86, (end-k));
#if defined(__SSE3__) && !defined(__AVX__)
			long long imask = 0x7FFFFFFFFFFFFFFF;
			double* fmask = (double*)&imask;

			for (size_t i = 0; i < source_size; i+=12)
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

			for (size_t i = 0; i < source_size; i+=24)
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
			for (size_t i = 0; i < source_size; i++)
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
#ifdef _OPENMP
	}
#endif

	return myTimer->stop();
}

double OperationMultipleEvalIterativeX86SimdModLinear::multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	size_t result_size = result.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();
    result.setAll(0.0);
    double* ptrAlpha = alpha.getPointer();
    double* ptrData = this->dataset_->getPointer();
    double* ptrResult = result.getPointer();
    double* ptrLevel = this->level_->getPointer();
    double* ptrIndex = this->index_->getPointer();

    if (this->dataset_->getNcols() % CHUNKDATAPOINTS_X86 != 0 || result_size != this->dataset_->getNcols())
    {
    	throw sg::base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

    myTimer->start();

#ifdef _OPENMP
    #pragma omp parallel
	{
		size_t chunksize = (result_size/omp_get_num_threads())+1;
		// assure that every subarray is 32-byte aligned
		if (chunksize % CHUNKDATAPOINTS_X86 != 0)
		{
			size_t remainder = chunksize % CHUNKDATAPOINTS_X86;
			size_t patch = CHUNKDATAPOINTS_X86 - remainder;
			chunksize += patch;
		}
    	size_t start = chunksize*omp_get_thread_num();
    	size_t end = std::min<size_t>(start+chunksize, result_size);
#else
    	size_t start = 0;
    	size_t end = result_size;
#endif
		for(size_t c = start; c < end; c+=std::min<size_t>((size_t)CHUNKDATAPOINTS_X86, (end-c)))
		{
			size_t data_end = std::min<size_t>((size_t)CHUNKDATAPOINTS_X86+c, end);

#ifdef __ICC
			#pragma ivdep
			#pragma vector aligned
#endif
			for (size_t i = c; i < data_end; i++)
			{
				ptrResult[i] = 0.0;
			}

			for (size_t m = 0; m < storageSize; m+=std::min<size_t>((size_t)CHUNKGRIDPOINTS_x86, (storageSize-m)))
			{
#if defined(__SSE3__) && !defined(__AVX__)
				size_t grid_inc = std::min<size_t>((size_t)CHUNKGRIDPOINTS_x86, (storageSize-m));

				long long imask = 0x7FFFFFFFFFFFFFFF;
				double* fmask = (double*)&imask;

				for (size_t i = c; i < c+CHUNKDATAPOINTS_X86; i+=12)
				{
					for (size_t j = m; j < m+grid_inc; j++)
					{
						__m128d support_0 = _mm_loaddup_pd(&(ptrAlpha[j]));
						__m128d support_1 = _mm_loaddup_pd(&(ptrAlpha[j]));
						__m128d support_2 = _mm_loaddup_pd(&(ptrAlpha[j]));
						__m128d support_3 = _mm_loaddup_pd(&(ptrAlpha[j]));
						__m128d support_4 = _mm_loaddup_pd(&(ptrAlpha[j]));
						__m128d support_5 = _mm_loaddup_pd(&(ptrAlpha[j]));

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
								__m128d eval_0 = _mm_load_pd(&(ptrData[(d*result_size)+i]));
								__m128d eval_1 = _mm_load_pd(&(ptrData[(d*result_size)+i+2]));
								__m128d eval_2 = _mm_load_pd(&(ptrData[(d*result_size)+i+4]));
								__m128d eval_3 = _mm_load_pd(&(ptrData[(d*result_size)+i+6]));
								__m128d eval_4 = _mm_load_pd(&(ptrData[(d*result_size)+i+8]));
								__m128d eval_5 = _mm_load_pd(&(ptrData[(d*result_size)+i+10]));

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
								__m128d eval_0 = _mm_load_pd(&(ptrData[(d*result_size)+i]));
								__m128d eval_1 = _mm_load_pd(&(ptrData[(d*result_size)+i+2]));
								__m128d eval_2 = _mm_load_pd(&(ptrData[(d*result_size)+i+4]));
								__m128d eval_3 = _mm_load_pd(&(ptrData[(d*result_size)+i+6]));
								__m128d eval_4 = _mm_load_pd(&(ptrData[(d*result_size)+i+8]));
								__m128d eval_5 = _mm_load_pd(&(ptrData[(d*result_size)+i+10]));

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
								__m128d eval_0 = _mm_load_pd(&(ptrData[(d*result_size)+i]));
								__m128d eval_1 = _mm_load_pd(&(ptrData[(d*result_size)+i+2]));
								__m128d eval_2 = _mm_load_pd(&(ptrData[(d*result_size)+i+4]));
								__m128d eval_3 = _mm_load_pd(&(ptrData[(d*result_size)+i+6]));
								__m128d eval_4 = _mm_load_pd(&(ptrData[(d*result_size)+i+8]));
								__m128d eval_5 = _mm_load_pd(&(ptrData[(d*result_size)+i+10]));

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
#endif
#if defined(__SSE3__) && defined(__AVX__)
				size_t grid_inc = std::min<size_t>((size_t)CHUNKGRIDPOINTS_x86, (storageSize-m));

				long long imask = 0x7FFFFFFFFFFFFFFF;
				double* fmask = (double*)&imask;

				for (size_t i = c; i < c+CHUNKDATAPOINTS_X86; i+=24)
				{
					for (size_t j = m; j < m+grid_inc; j++)
					{
						__m256d support_0 = _mm256_broadcast_sd(&(ptrAlpha[j]));
						__m256d support_1 = _mm256_broadcast_sd(&(ptrAlpha[j]));
						__m256d support_2 = _mm256_broadcast_sd(&(ptrAlpha[j]));
						__m256d support_3 = _mm256_broadcast_sd(&(ptrAlpha[j]));
						__m256d support_4 = _mm256_broadcast_sd(&(ptrAlpha[j]));
						__m256d support_5 = _mm256_broadcast_sd(&(ptrAlpha[j]));

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
								__m256d eval_0 = _mm256_load_pd(&(ptrData[(d*result_size)+i]));
								__m256d eval_1 = _mm256_load_pd(&(ptrData[(d*result_size)+i+4]));
								__m256d eval_2 = _mm256_load_pd(&(ptrData[(d*result_size)+i+8]));
								__m256d eval_3 = _mm256_load_pd(&(ptrData[(d*result_size)+i+12]));
								__m256d eval_4 = _mm256_load_pd(&(ptrData[(d*result_size)+i+16]));
								__m256d eval_5 = _mm256_load_pd(&(ptrData[(d*result_size)+i+20]));

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
							else if (ptrIndex[(j*dims)+d] == (ptrLevel[(j*dims)+d] - 1.0))
							{
								__m256d eval_0 = _mm256_load_pd(&(ptrData[(d*result_size)+i]));
								__m256d eval_1 = _mm256_load_pd(&(ptrData[(d*result_size)+i+4]));
								__m256d eval_2 = _mm256_load_pd(&(ptrData[(d*result_size)+i+8]));
								__m256d eval_3 = _mm256_load_pd(&(ptrData[(d*result_size)+i+12]));
								__m256d eval_4 = _mm256_load_pd(&(ptrData[(d*result_size)+i+16]));
								__m256d eval_5 = _mm256_load_pd(&(ptrData[(d*result_size)+i+20]));

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
								__m256d eval_0 = _mm256_load_pd(&(ptrData[(d*result_size)+i]));
								__m256d eval_1 = _mm256_load_pd(&(ptrData[(d*result_size)+i+4]));
								__m256d eval_2 = _mm256_load_pd(&(ptrData[(d*result_size)+i+8]));
								__m256d eval_3 = _mm256_load_pd(&(ptrData[(d*result_size)+i+12]));
								__m256d eval_4 = _mm256_load_pd(&(ptrData[(d*result_size)+i+16]));
								__m256d eval_5 = _mm256_load_pd(&(ptrData[(d*result_size)+i+20]));

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
#endif
#if !defined(__SSE3__) && !defined(__AVX__)
				size_t grid_end = std::min<size_t>((size_t)CHUNKGRIDPOINTS_x86+m, storageSize);

				for (size_t i = c; i < data_end; i++)
				{
					for (size_t j = m; j < grid_end; j++)
					{
						double curSupport = ptrAlpha[j];

						for (size_t d = 0; d < dims; d++)
						{
							if (ptrLevel[(j*dims)+d] == 2.0)
							{
								// nothing to do (mult with 1)
							}
							else if (ptrIndex[(j*dims)+d] == 1.0)
							{
								double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*result_size)+i]));
								eval = 2.0 - eval;
								double localSupport = std::max<double>(eval, 0.0);
								curSupport *= localSupport;
							}
							else if (ptrIndex[(j*dims)+d] == (ptrLevel[(j*dims)+d] - 1.0))
							{
								double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*result_size)+i]));
								double index_calc = eval - (ptrIndex[(j*dims)+d]);
								double last = 1.0 + index_calc;
								double localSupport = std::max<double>(last, 0.0);
								curSupport *= localSupport;
							}
							else
							{
								double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*result_size)+i]));
								double index_calc = eval - (ptrIndex[(j*dims)+d]);
								double abs = fabs(index_calc);
								double last = 1.0 - abs;
								double localSupport = std::max<float>(last, 0.0);
								curSupport *= localSupport;
							}
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

}
