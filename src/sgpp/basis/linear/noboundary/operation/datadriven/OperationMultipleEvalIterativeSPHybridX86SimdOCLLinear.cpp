/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear.hpp"
#include "exception/operation_exception.hpp"
#include "tools/common/AlignedMemory.hpp"

#ifdef __ICC
#include <x86intrin.h>
#endif

namespace sg
{
namespace parallel
{

OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear::OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear(sg::base::GridStorage* storage, sg::base::DataMatrixSP* dataset) : sg::base::OperationMultipleEvalVectorizedSP(dataset)
{
	this->storage = storage;

	this->level_ = new sg::base::DataMatrixSP(storage->size(), storage->dim());
	this->index_ = new sg::base::DataMatrixSP(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myTimer = new sg::base::SGppStopwatch();
	myOCLKernels = new OCLKernels();

	_tuningMult = new sg::base::TwoPartitionAutoTuning(dataset->getNrows(), 128);
	_tuningMultTrans = new sg::base::TwoPartitionAutoTuning(storage->size(), 128);
}

OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear::~OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear()
{
	delete myTimer;
	delete myOCLKernels;
	delete _tuningMult;
	delete _tuningMultTrans;
}

void OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear::rebuildLevelAndIndex()
{
	delete this->level_;
	delete this->index_;

	this->level_ = new sg::base::DataMatrixSP(storage->size(), storage->dim());
	this->index_ = new sg::base::DataMatrixSP(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myOCLKernels->resetKernels();

	_tuningMultTrans->setProblemSize(storage->size());
}

double OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear::multTransposeVectorized(sg::base::DataVectorSP& source, sg::base::DataVectorSP& result)
{
	size_t source_size = source.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    double gpu_time = 0.0;
    double cpu_time = 0.0;

    result.setAll(0.0f);

    float* ptrSource = source.getPointer();
    float* ptrData = this->dataset_->getPointer();
    float* ptrLevel = this->level_->getPointer();
    float* ptrIndex = this->index_->getPointer();
    float* ptrGlobalResult = result.getPointer();

    if (this->dataset_->getNrows() % 128 != 0 || source_size != this->dataset_->getNrows())
    {
    	throw sg::base::operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
    }

    // split result into GPU and CPU partition
    size_t gpu_partition = storageSize - _tuningMultTrans->getPartition1Size();

    // Do on-demand transpose
	float* ptrTransData = new float[dims*source_size];

	#pragma omp parallel shared(gpu_time, cpu_time)
    {
		#pragma omp single nowait
    	{
			#pragma omp task shared(gpu_time, cpu_time)
    		{
    			if (gpu_partition > 0)
    			{
    				gpu_time = myOCLKernels->multTransSPOCL(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims, gpu_partition);
    			}
    		}

			#pragma omp task shared(gpu_time, cpu_time)
    		{
    			myTimer->start();
#ifdef __ICC
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
						__m128 res = _mm_set1_ps(0.0f);
						int imask = 0x7FFFFFFF;
						float* fmask = (float*)&imask;

						for (size_t i = 0; i < source_size; i+=16)
						{
							__m128 support_0 = _mm_load_ps(&(ptrSource[i+0]));
							__m128 support_1 = _mm_load_ps(&(ptrSource[i+4]));
							__m128 support_2 = _mm_load_ps(&(ptrSource[i+8]));
							__m128 support_3 = _mm_load_ps(&(ptrSource[i+12]));

							__m128 one = _mm_set1_ps(1.0f);
							__m128 zero = _mm_set1_ps(0.0f);
							__m128 mask = _mm_set1_ps(*fmask);

							for (size_t d = 0; d < dims; d++)
							{
								__m128 eval_0 = _mm_load_ps(&(ptrTransData[(d*source_size)+i]));
								__m128 eval_1 = _mm_load_ps(&(ptrTransData[(d*source_size)+i+4]));
								__m128 eval_2 = _mm_load_ps(&(ptrTransData[(d*source_size)+i+8]));
								__m128 eval_3 = _mm_load_ps(&(ptrTransData[(d*source_size)+i+12]));;

								__m128 level = _mm_load1_ps(&(ptrLevel[(j*dims)+d]));
								__m128 index = _mm_load1_ps(&(ptrIndex[(j*dims)+d]));

								eval_0 = _mm_mul_ps(eval_0, level);
								eval_1 = _mm_mul_ps(eval_1, level);
								eval_2 = _mm_mul_ps(eval_2, level);
								eval_3 = _mm_mul_ps(eval_3, level);

								eval_0 = _mm_sub_ps(eval_0, index);
								eval_1 = _mm_sub_ps(eval_1, index);
								eval_2 = _mm_sub_ps(eval_2, index);
								eval_3 = _mm_sub_ps(eval_3, index);

								eval_0 = _mm_and_ps(mask, eval_0);
								eval_1 = _mm_and_ps(mask, eval_1);
								eval_2 = _mm_and_ps(mask, eval_2);
								eval_3 = _mm_and_ps(mask, eval_3);

								eval_0 = _mm_sub_ps(one, eval_0);
								eval_1 = _mm_sub_ps(one, eval_1);
								eval_2 = _mm_sub_ps(one, eval_2);
								eval_3 = _mm_sub_ps(one, eval_3);

								eval_0 = _mm_max_ps(zero, eval_0);
								eval_1 = _mm_max_ps(zero, eval_1);
								eval_2 = _mm_max_ps(zero, eval_2);
								eval_3 = _mm_max_ps(zero, eval_3);

								support_0 = _mm_mul_ps(support_0, eval_0);
								support_1 = _mm_mul_ps(support_1, eval_1);
								support_2 = _mm_mul_ps(support_2, eval_2);
								support_3 = _mm_mul_ps(support_3, eval_3);
							}

							support_0 = _mm_add_ps(support_0, support_1);
							support_2 = _mm_add_ps(support_2, support_3);
							support_0 = _mm_add_ps(support_0, support_2);

							res = _mm_add_ps(res, support_0);
						}

						res = _mm_hadd_ps(res, res);
						res = _mm_hadd_ps(res, res);

						_mm_store_ss(&(ptrGlobalResult[j]), res);
					}
				}
#else
				for (size_t j = gpu_partition; j < storageSize; j++)
				{
					#pragma omp task firstprivate(j)
					{
						ptrGlobalResult[j] = 0.0f;

						for (size_t i = 0; i < source_size; i++)
						{
							float curSupport = ptrSource[i];

							for (size_t d = 0; d < dims; d++)
							{
								float eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
								float index_calc = eval - (ptrIndex[(j*dims)+d]);
								float abs = fabs(index_calc);
								float last = 1.0f - abs;
								float localSupport = std::max<float>(last, 0.0f);
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

    _tuningMultTrans->setExecutionTimes(cpu_time, gpu_time);

    double time = 0.0;
    //cleanup
    delete[] ptrTransData;

	return time;
}

double OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear::multVectorized(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result)
{
	size_t result_size = result.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    double gpu_time = 0.0;
    double cpu_time = 0.0;

    result.setAll(0.0f);

    float* ptrAlpha = alpha.getPointer();
    float* ptrData = this->dataset_->getPointer();
    float* ptrResult = result.getPointer();
    float* ptrLevel = this->level_->getPointer();
    float* ptrIndex = this->index_->getPointer();

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
    				gpu_time = myOCLKernels->multSPOCL(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims, gpu_partition);
    			}
    		}

			#pragma omp task shared(gpu_time, cpu_time)
    		{
    			myTimer->start();
#ifdef __ICC
				for (size_t i = gpu_partition; i < result_size; i+=16)
				{
					#pragma omp task firstprivate(i)
					{
						int imask = 0x7FFFFFFF;
						float* fmask = (float*)&imask;

						__m128 res_0 = _mm_load_ps(&(ptrResult[i]));
						__m128 res_1 = _mm_load_ps(&(ptrResult[i+4]));
						__m128 res_2 = _mm_load_ps(&(ptrResult[i+8]));
						__m128 res_3 = _mm_load_ps(&(ptrResult[i+12]));

						// Do on-demand transpose
						float* ptrTransData = new float[dims*16];
						for (size_t n = 0; n < 16; n++)
						{
							for(size_t d = 0; d < dims; d++)
							{
								ptrTransData[(d*16)+n] = ptrData[((i+n)*dims)+d];
							}
						}

						for (size_t j = 0; j < storageSize; j++)
						{
							__m128 support_0 = _mm_load1_ps(&(ptrAlpha[j]));
							__m128 support_1 = _mm_load1_ps(&(ptrAlpha[j]));
							__m128 support_2 = _mm_load1_ps(&(ptrAlpha[j]));
							__m128 support_3 = _mm_load1_ps(&(ptrAlpha[j]));

							__m128 one = _mm_set1_ps(1.0f);
							__m128 zero = _mm_set1_ps(0.0f);
							__m128 mask = _mm_set1_ps(*fmask);

							for (size_t d = 0; d < dims; d++)
							{
								__m128 eval_0 = _mm_load_ps(&(ptrTransData[(d*16)]));
								__m128 eval_1 = _mm_load_ps(&(ptrTransData[(d*16)+4]));
								__m128 eval_2 = _mm_load_ps(&(ptrTransData[(d*16)+8]));
								__m128 eval_3 = _mm_load_ps(&(ptrTransData[(d*16)+12]));;

								__m128 level = _mm_load1_ps(&(ptrLevel[(j*dims)+d]));
								__m128 index = _mm_load1_ps(&(ptrIndex[(j*dims)+d]));

								eval_0 = _mm_mul_ps(eval_0, level);
								eval_1 = _mm_mul_ps(eval_1, level);
								eval_2 = _mm_mul_ps(eval_2, level);
								eval_3 = _mm_mul_ps(eval_3, level);

								eval_0 = _mm_sub_ps(eval_0, index);
								eval_1 = _mm_sub_ps(eval_1, index);
								eval_2 = _mm_sub_ps(eval_2, index);
								eval_3 = _mm_sub_ps(eval_3, index);

								eval_0 = _mm_and_ps(mask, eval_0);
								eval_1 = _mm_and_ps(mask, eval_1);
								eval_2 = _mm_and_ps(mask, eval_2);
								eval_3 = _mm_and_ps(mask, eval_3);

								eval_0 = _mm_sub_ps(one, eval_0);
								eval_1 = _mm_sub_ps(one, eval_1);
								eval_2 = _mm_sub_ps(one, eval_2);
								eval_3 = _mm_sub_ps(one, eval_3);

								eval_0 = _mm_max_ps(zero, eval_0);
								eval_1 = _mm_max_ps(zero, eval_1);
								eval_2 = _mm_max_ps(zero, eval_2);
								eval_3 = _mm_max_ps(zero, eval_3);

								support_0 = _mm_mul_ps(support_0, eval_0);
								support_1 = _mm_mul_ps(support_1, eval_1);
								support_2 = _mm_mul_ps(support_2, eval_2);
								support_3 = _mm_mul_ps(support_3, eval_3);
							}

							res_0 = _mm_add_ps(res_0, support_0);
							res_1 = _mm_add_ps(res_1, support_1);
							res_2 = _mm_add_ps(res_2, support_2);
							res_3 = _mm_add_ps(res_3, support_3);
						}

						delete[] ptrTransData;

						_mm_store_ps(&(ptrResult[i]), res_0);
						_mm_store_ps(&(ptrResult[i+4]), res_1);
						_mm_store_ps(&(ptrResult[i+8]), res_2);
						_mm_store_ps(&(ptrResult[i+12]), res_3);
					}
				}
#else
				for (size_t i = gpu_partition; i < result_size; i++)
				{
					#pragma omp task firstprivate(i)
					{
						for (size_t j = 0; j < storageSize; j++)
						{
							float curSupport = ptrAlpha[j];

							for (size_t d = 0; d < dims; d++)
							{
								float eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
								float index_calc = eval - (ptrIndex[(j*dims)+d]);
								float abs = fabs(index_calc);
								float last = 1.0f - abs;
								float localSupport = std::max<double>(last, 0.0f);
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

    std::cout << cpu_time << " " << gpu_time << std::endl;

    _tuningMult->setExecutionTimes(cpu_time, gpu_time);

    double time = 0.0;

   	return time;
}

}
}
