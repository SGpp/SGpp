/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeHybridX86SimdOCLLinear.hpp"
#include "base/exception/operation_exception.hpp"
#include "base/tools/AlignedMemory.hpp"

#if defined(__SSE3__) || defined(__AVX__)
#include <x86intrin.h>
#endif

#ifdef __USEAVX128__
#undef __AVX__
#endif

namespace sg
{
namespace parallel
{

OperationMultipleEvalIterativeHybridX86SimdOCLLinear::OperationMultipleEvalIterativeHybridX86SimdOCLLinear(sg::base::GridStorage* storage, sg::base::DataMatrix* dataset) : sg::base::OperationMultipleEvalVectorized(dataset)
{
	this->storage = storage;

	this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
	this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myTimer = new sg::base::SGppStopwatch();
	myOCLKernels = new OCLKernels();

	_tuningMult = new sg::base::TwoPartitionAutoTuning(dataset->getNrows(), 128, 10, 0.75, 3);
	_tuningMultTrans = new sg::base::TwoPartitionAutoTuning(storage->size(), 128, 10, 0.85, 3);
}

OperationMultipleEvalIterativeHybridX86SimdOCLLinear::~OperationMultipleEvalIterativeHybridX86SimdOCLLinear()
{
	delete myTimer;
	delete myOCLKernels;
	delete _tuningMult;
	delete _tuningMultTrans;
}

void OperationMultipleEvalIterativeHybridX86SimdOCLLinear::rebuildLevelAndIndex()
{
	delete this->level_;
	delete this->index_;

	this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
	this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myOCLKernels->resetKernels();

	_tuningMultTrans->setProblemSize(storage->size());
	_tuningMult->softResetAutoTuning();
}

double OperationMultipleEvalIterativeHybridX86SimdOCLLinear::multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result)
{
	size_t source_size = source.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    double gpu_time = 0.0;
    double cpu_time = 0.0;

    result.setAll(0.0);

    double* ptrSource = source.getPointer();
    double* ptrData = this->dataset_->getPointer();
    double* ptrLevel = this->level_->getPointer();
    double* ptrIndex = this->index_->getPointer();
    double* ptrGlobalResult = result.getPointer();

    if (this->dataset_->getNrows() % 128 != 0 || source_size != this->dataset_->getNrows())
    {
    	throw sg::base::operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
    }

    // split result into GPU and CPU partition
    size_t gpu_partition = storageSize - _tuningMultTrans->getPartition1Size();

    // Do on-demand transpose
	double* ptrTransData = new double[dims*source_size];

	#pragma omp parallel shared(gpu_time, cpu_time)
    {
		#pragma omp single nowait
    	{
			#pragma omp task shared(gpu_time, cpu_time)
    		{
    			if (gpu_partition > 0)
    			{
    				gpu_time = myOCLKernels->multTransOCL(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims, gpu_partition);
    			}
    		}

			#pragma omp task shared(gpu_time, cpu_time)
    		{
    			myTimer->start();
#if defined(__SSE3__) && !defined(__AVX__)
    			for (size_t n = 0; n < source_size; n++)
    			{
    				for(size_t d = 0; d < dims; d++)
    				{
    					ptrTransData[(d*source_size)+n] = ptrData[(n*dims)+d];
    				}
    			}

    			for (size_t j = gpu_partition; j < storageSize; j++)
				{
					#pragma omp task firstprivate(j)
					{
						long long imask = 0x7FFFFFFFFFFFFFFF;
						double* fmask = (double*)&imask;

						__m128d res = _mm_set1_pd(0.0f);

						for (size_t i = 0; i < source_size; i+=8)
						{
							__m128d support_0 = _mm_load_pd(&(ptrSource[i+0]));
							__m128d support_1 = _mm_load_pd(&(ptrSource[i+2]));
							__m128d support_2 = _mm_load_pd(&(ptrSource[i+4]));
							__m128d support_3 = _mm_load_pd(&(ptrSource[i+6]));

							__m128d one = _mm_set1_pd(1.0);
							__m128d zero = _mm_set1_pd(0.0);
							__m128d mask = _mm_set1_pd(*fmask);

							for (size_t d = 0; d < dims; d++)
							{
								__m128d eval_0 = _mm_load_pd(&(ptrTransData[(d*source_size)+i+0]));
								__m128d eval_1 = _mm_load_pd(&(ptrTransData[(d*source_size)+i+2]));
								__m128d eval_2 = _mm_load_pd(&(ptrTransData[(d*source_size)+i+4]));
								__m128d eval_3 = _mm_load_pd(&(ptrTransData[(d*source_size)+i+6]));;

								__m128d level = _mm_loaddup_pd(&(ptrLevel[(j*dims)+d]));
								__m128d index = _mm_loaddup_pd(&(ptrIndex[(j*dims)+d]));

#ifdef __FMA4__
								eval_0 = _mm_msub_pd(eval_0, level, index);
								eval_1 = _mm_msub_pd(eval_1, level, index);
								eval_2 = _mm_msub_pd(eval_2, level, index);
								eval_3 = _mm_msub_pd(eval_3, level, index);
#else
								eval_0 = _mm_sub_pd(_mm_mul_pd(eval_0, level), index);
								eval_1 = _mm_sub_pd(_mm_mul_pd(eval_1, level), index);
								eval_2 = _mm_sub_pd(_mm_mul_pd(eval_2, level), index);
								eval_3 = _mm_sub_pd(_mm_mul_pd(eval_3, level), index);
#endif

								eval_0 = _mm_and_pd(mask, eval_0);
								eval_1 = _mm_and_pd(mask, eval_1);
								eval_2 = _mm_and_pd(mask, eval_2);
								eval_3 = _mm_and_pd(mask, eval_3);

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

							support_0 = _mm_add_pd(support_0, support_1);
							support_2 = _mm_add_pd(support_2, support_3);
							support_0 = _mm_add_pd(support_0, support_2);

							res = _mm_add_pd(res, support_0);
						}

						res = _mm_hadd_pd(res, res);

						_mm_store_sd(&(ptrGlobalResult[j]), res);
					}
				}
#endif
#if defined(__SSE3__) && defined(__AVX__)
    			for (size_t n = 0; n < source_size; n++)
    			{
    				for(size_t d = 0; d < dims; d++)
    				{
    					ptrTransData[(d*source_size)+n] = ptrData[(n*dims)+d];
    				}
    			}

    			for (size_t j = gpu_partition; j < storageSize; j++)
				{
					#pragma omp task firstprivate(j)
					{
						long long imask = 0x7FFFFFFFFFFFFFFF;
						double* fmask = (double*)&imask;

						__m256d res = _mm_set1_pd(0.0f);

						for (size_t i = 0; i < source_size; i+=16)
						{
							__m256d support_0 = _mm256_load_pd(&(ptrSource[i+0]));
							__m256d support_1 = _mm256_load_pd(&(ptrSource[i+4]));
							__m256d support_2 = _mm256_load_pd(&(ptrSource[i+8]));
							__m256d support_3 = _mm256_load_pd(&(ptrSource[i+12]));

							__m256d one = _mm256_set1_pd(1.0);
							__m256d zero = _mm256_set1_pd(0.0);
							__m256d mask = _mm256_set1_pd(*fmask);

							for (size_t d = 0; d < dims; d++)
							{
								__m256d eval_0 = _mm256_load_pd(&(ptrTransData[(d*source_size)+i+0]));
								__m256d eval_1 = _mm256_load_pd(&(ptrTransData[(d*source_size)+i+4]));
								__m256d eval_2 = _mm256_load_pd(&(ptrTransData[(d*source_size)+i+8]));
								__m256d eval_3 = _mm256_load_pd(&(ptrTransData[(d*source_size)+i+12]));;

								__m256d level = _mm256_broadcast_sd(&(ptrLevel[(j*dims)+d]));
								__m256d index = _mm256_broadcast_sd(&(ptrIndex[(j*dims)+d]));
#ifdef __FMA4__
								eval_0 = _mm256_msub_pd(eval_0, level, index);
								eval_1 = _mm256_msub_pd(eval_1, level, index);
								eval_2 = _mm256_msub_pd(eval_2, level, index);
								eval_3 = _mm256_msub_pd(eval_3, level, index);
#else
								eval_0 = _mm256_sub_pd(_mm256_mul_pd(eval_0, level), index);
								eval_1 = _mm256_sub_pd(_mm256_mul_pd(eval_1, level), index);
								eval_2 = _mm256_sub_pd(_mm256_mul_pd(eval_2, level), index);
								eval_3 = _mm256_sub_pd(_mm256_mul_pd(eval_3, level), index);
#endif
								eval_0 = _mm256_and_pd(mask, eval_0);
								eval_1 = _mm256_and_pd(mask, eval_1);
								eval_2 = _mm256_and_pd(mask, eval_2);
								eval_3 = _mm256_and_pd(mask, eval_3);

								eval_0 = _mm256_sub_pd(one, eval_0);
								eval_1 = _mm256_sub_pd(one, eval_1);
								eval_2 = _mm256_sub_pd(one, eval_2);
								eval_3 = _mm256_sub_pd(one, eval_3);

								eval_0 = _mm256_max_pd(zero, eval_0);
								eval_1 = _mm256_max_pd(zero, eval_1);
								eval_2 = _mm256_max_pd(zero, eval_2);
								eval_3 = _mm256_max_pd(zero, eval_3);

								support_0 = _mm256_mul_pd(support_0, eval_0);
								support_1 = _mm256_mul_pd(support_1, eval_1);
								support_2 = _mm256_mul_pd(support_2, eval_2);
								support_3 = _mm256_mul_pd(support_3, eval_3);
							}

							support_0 = _mm256_add_pd(support_0, support_1);
							support_2 = _mm256_add_pd(support_2, support_3);
							support_0 = _mm256_add_pd(support_0, support_2);

							res = _mm256_add_pd(res, support_0);
						}

						const __m256i ldStMaskAVX = _mm256_set_epi64x(0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xFFFFFFFFFFFFFFFF);

						res = _mm256_hadd_pd(res, res);
						__m256d tmp = _mm256_permute2f128_pd(res, res, 0x81);
						res = _mm256_add_pd(res, tmp);

						_mm256_maskstore_pd(&(ptrResult[j]), ldStMaskAVX, res_0);
					}
				}
#endif
#if !defined(__SSE3__) && !defined(__AVX__)
				for (size_t j = gpu_partition; j < storageSize; j++)
				{
					#pragma omp task firstprivate(j)
					{
						ptrGlobalResult[j] = 0.0;

						for (size_t i = 0; i < source_size; i++)
						{
							double curSupport = ptrSource[i];

							for (size_t d = 0; d < dims; d++)
							{
								double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
								double index_calc = eval - (ptrIndex[(j*dims)+d]);
								double abs = fabs(index_calc);
								double last = 1.0 - abs;
								double localSupport = std::max<double>(last, 0.0);
								curSupport *= localSupport;
							}

							ptrGlobalResult[j] += curSupport;
						}
					}
				}
#endif
				cpu_time = myTimer->stop();
    		}

			#pragma omp taskwait
    	}
    }

    double time = 0.0;

    _tuningMultTrans->setExecutionTimes(cpu_time, gpu_time);

    //cleanup
    delete[] ptrTransData;

	return time;
}

double OperationMultipleEvalIterativeHybridX86SimdOCLLinear::multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	size_t result_size = result.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    double gpu_time = 0.0;
    double cpu_time = 0.0;

    result.setAll(0.0);

    double* ptrAlpha = alpha.getPointer();
    double* ptrData = this->dataset_->getPointer();
    double* ptrResult = result.getPointer();
    double* ptrLevel = this->level_->getPointer();
    double* ptrIndex = this->index_->getPointer();

    if (this->dataset_->getNrows() % 128 != 0 || result_size != this->dataset_->getNrows())
    {
    	throw sg::base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

    // split result into GPU and CPU partition
    size_t gpu_partition = result_size - _tuningMult->getPartition1Size();

	#pragma omp parallel shared(gpu_time, cpu_time)
    {
		#pragma omp single nowait
    	{
			#pragma omp task shared(gpu_time, cpu_time)
    		{
    			if (gpu_partition > 0)
    			{
    				gpu_time = myOCLKernels->multOCL(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims, gpu_partition);
    			}
    		}

			#pragma omp task shared(gpu_time, cpu_time)
    		{
    			myTimer->start();
#if defined(__SSE3__) && !defined(__AVX__)
				for (size_t i = gpu_partition; i < result_size; i+=8)
				{
					#pragma omp task firstprivate(i)
					{
						long long imask = 0x7FFFFFFFFFFFFFFF;
						double* fmask = (double*)&imask;

						__m128d res_0 = _mm_load_pd(&(ptrResult[i+0]));
						__m128d res_1 = _mm_load_pd(&(ptrResult[i+2]));
						__m128d res_2 = _mm_load_pd(&(ptrResult[i+4]));
						__m128d res_3 = _mm_load_pd(&(ptrResult[i+6]));

						// Do on-demand transpose
						double* ptrTransData = new double[dims*8];
						for (size_t n = 0; n < 8; n++)
						{
							for(size_t d = 0; d < dims; d++)
							{
								ptrTransData[(d*8)+n] = ptrData[((i+n)*dims)+d];
							}
						}

						for (size_t j = 0; j < storageSize; j++)
						{
							__m128d support_0 = _mm_loaddup_pd(&(ptrAlpha[j]));
							__m128d support_1 = _mm_loaddup_pd(&(ptrAlpha[j]));
							__m128d support_2 = _mm_loaddup_pd(&(ptrAlpha[j]));
							__m128d support_3 = _mm_loaddup_pd(&(ptrAlpha[j]));

							__m128d one = _mm_set1_pd(1.0);
							__m128d zero = _mm_set1_pd(0.0);
							__m128d mask = _mm_set1_pd(*fmask);

							for (size_t d = 0; d < dims; d++)
							{
								__m128d eval_0 = _mm_load_pd(&(ptrTransData[(d*8)+0]));
								__m128d eval_1 = _mm_load_pd(&(ptrTransData[(d*8)+2]));
								__m128d eval_2 = _mm_load_pd(&(ptrTransData[(d*8)+4]));
								__m128d eval_3 = _mm_load_pd(&(ptrTransData[(d*8)+6]));;

								__m128d level = _mm_loaddup_pd(&(ptrLevel[(j*dims)+d]));
								__m128d index = _mm_loaddup_pd(&(ptrIndex[(j*dims)+d]));
#ifdef __FMA4__
								eval_0 = _mm_msub_pd(eval_0, level, index);
								eval_1 = _mm_msub_pd(eval_1, level, index);
								eval_2 = _mm_msub_pd(eval_2, level, index);
								eval_3 = _mm_msub_pd(eval_3, level, index);
#else
								eval_0 = _mm_sub_pd(_mm_mul_pd(eval_0, level), index);
								eval_1 = _mm_sub_pd(_mm_mul_pd(eval_1, level), index);
								eval_2 = _mm_sub_pd(_mm_mul_pd(eval_2, level), index);
								eval_3 = _mm_sub_pd(_mm_mul_pd(eval_3, level), index);
#endif
								eval_0 = _mm_and_pd(mask, eval_0);
								eval_1 = _mm_and_pd(mask, eval_1);
								eval_2 = _mm_and_pd(mask, eval_2);
								eval_3 = _mm_and_pd(mask, eval_3);

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

							res_0 = _mm_add_pd(res_0, support_0);
							res_1 = _mm_add_pd(res_1, support_1);
							res_2 = _mm_add_pd(res_2, support_2);
							res_3 = _mm_add_pd(res_3, support_3);
						}

						delete[] ptrTransData;

						_mm_store_pd(&(ptrResult[i+0]), res_0);
						_mm_store_pd(&(ptrResult[i+2]), res_1);
						_mm_store_pd(&(ptrResult[i+4]), res_2);
						_mm_store_pd(&(ptrResult[i+6]), res_3);
					}
				}
#endif
#if defined(__SSE3__) && defined(__AVX__)
				for (size_t i = gpu_partition; i < result_size; i+=16)
				{
					#pragma omp task firstprivate(i)
					{
						long long imask = 0x7FFFFFFFFFFFFFFF;
						double* fmask = (double*)&imask;

						__m256d res_0 = _mm256_load_pd(&(ptrResult[i+0]));
						__m256d res_1 = _mm256_load_pd(&(ptrResult[i+4]));
						__m256d res_2 = _mm256_load_pd(&(ptrResult[i+8]));
						__m256d res_3 = _mm256_load_pd(&(ptrResult[i+12]));

						// Do on-demand transpose
						double* ptrTransData = new double[dims*16];
						for (size_t n = 0; n < 16; n++)
						{
							for(size_t d = 0; d < dims; d++)
							{
								ptrTransData[(d*16)+n] = ptrData[((i+n)*dims)+d];
							}
						}

						for (size_t j = 0; j < storageSize; j++)
						{
							__m256d support_0 = _mm256_broadcast_sd(&(ptrAlpha[j]));
							__m256d support_1 = _mm256_broadcast_sd(&(ptrAlpha[j]));
							__m256d support_2 = _mm256_broadcast_sd(&(ptrAlpha[j]));
							__m256d support_3 = _mm256_broadcast_sd(&(ptrAlpha[j]));

							__m256d one = _mm256_set1_pd(1.0);
							__m256d zero = _mm256_set1_pd(0.0);
							__m256d mask = _mm256_set1_pd(*fmask);

							for (size_t d = 0; d < dims; d++)
							{
								__m256d eval_0 = _mm256_load_pd(&(ptrTransData[(d*16)+0]));
								__m256d eval_1 = _mm256_load_pd(&(ptrTransData[(d*16)+4]));
								__m256d eval_2 = _mm256_load_pd(&(ptrTransData[(d*16)+8]));
								__m256d eval_3 = _mm256_load_pd(&(ptrTransData[(d*16)+12]));;

								__m256d level = _mm256_broadcast_sd(&(ptrLevel[(j*dims)+d]));
								__m256d index = _mm256_broadcast_sd(&(ptrIndex[(j*dims)+d]));
#ifdef __FMA4__
								eval_0 = _mm256_msub_pd(eval_0, level, index);
								eval_1 = _mm256_msub_pd(eval_1, level, index);
								eval_2 = _mm256_msub_pd(eval_2, level, index);
								eval_3 = _mm256_msub_pd(eval_3, level, index);
#else
								eval_0 = _mm256_sub_pd(_mm256_mul_pd(eval_0, level), index);
								eval_1 = _mm256_sub_pd(_mm256_mul_pd(eval_1, level), index);
								eval_2 = _mm256_sub_pd(_mm256_mul_pd(eval_2, level), index);
								eval_3 = _mm256_sub_pd(_mm256_mul_pd(eval_3, level), index);
#endif
								eval_0 = _mm256_and_pd(mask, eval_0);
								eval_1 = _mm256_and_pd(mask, eval_1);
								eval_2 = _mm256_and_pd(mask, eval_2);
								eval_3 = _mm256_and_pd(mask, eval_3);

								eval_0 = _mm256_sub_pd(one, eval_0);
								eval_1 = _mm256_sub_pd(one, eval_1);
								eval_2 = _mm256_sub_pd(one, eval_2);
								eval_3 = _mm256_sub_pd(one, eval_3);

								eval_0 = _mm256_max_pd(zero, eval_0);
								eval_1 = _mm256_max_pd(zero, eval_1);
								eval_2 = _mm256_max_pd(zero, eval_2);
								eval_3 = _mm256_max_pd(zero, eval_3);

								support_0 = _mm256_mul_pd(support_0, eval_0);
								support_1 = _mm256_mul_pd(support_1, eval_1);
								support_2 = _mm256_mul_pd(support_2, eval_2);
								support_3 = _mm256_mul_pd(support_3, eval_3);
							}

							res_0 = _mm256_add_pd(res_0, support_0);
							res_1 = _mm256_add_pd(res_1, support_1);
							res_2 = _mm256_add_pd(res_2, support_2);
							res_3 = _mm256_add_pd(res_3, support_3);
						}

						delete[] ptrTransData;

						_mm256_store_pd(&(ptrResult[i+0]), res_0);
						_mm256_store_pd(&(ptrResult[i+4]), res_1);
						_mm256_store_pd(&(ptrResult[i+8]), res_2);
						_mm256_store_pd(&(ptrResult[i+12]), res_3);
					}
				}
#endif
#if !defined(__SSE3__) && !defined(__AVX__)
				for (size_t i = gpu_partition; i < result_size; i++)
				{
					#pragma omp task firstprivate(i)
					{
						for (size_t j = 0; j < storageSize; j++)
						{
							double curSupport = ptrAlpha[j];

							for (size_t d = 0; d < dims; d++)
							{
								double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
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
#endif
				cpu_time = myTimer->stop();
    		}

			#pragma omp taskwait
    	}
    }
    double time = 0.0;

    _tuningMult->setExecutionTimes(cpu_time, gpu_time);

   	return time;
}

}
}
