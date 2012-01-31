/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeSPOCLModLinear.hpp"
#include "base/exception/operation_exception.hpp"

namespace sg
{

namespace parallel
{

  OperationMultipleEvalIterativeSPOCLModLinear::OperationMultipleEvalIterativeSPOCLModLinear(base::GridStorage* storage, base::DataMatrixSP* dataset) : sg::parallel::OperationMultipleEvalVectorizedSP(dataset)
  {
    this->storage = storage;

    this->level_ = new base::DataMatrixSP(storage->size(), storage->dim());
    this->index_ = new base::DataMatrixSP(storage->size(), storage->dim());

    storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

    myTimer = new base::SGppStopwatch();
    myOCLKernels = new OCLKernels();
  }

  OperationMultipleEvalIterativeSPOCLModLinear::~OperationMultipleEvalIterativeSPOCLModLinear()
  {
    delete myTimer;
    delete myOCLKernels;
  }

  void OperationMultipleEvalIterativeSPOCLModLinear::rebuildLevelAndIndex()
  {
    delete this->level_;
    delete this->index_;

    this->level_ = new base::DataMatrixSP(storage->size(), storage->dim());
    this->index_ = new base::DataMatrixSP(storage->size(), storage->dim());

    storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

    myOCLKernels->resetKernels();
  }

  double OperationMultipleEvalIterativeSPOCLModLinear::multTransposeVectorized(base::DataVectorSP& source, base::DataVectorSP& result)
  {
    size_t source_size = source.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    result.setAll(0.0);

    float* ptrSource = source.getPointer();
    float* ptrData = this->dataset_->getPointer();
    float* ptrLevel = this->level_->getPointer();
    float* ptrIndex = this->index_->getPointer();
    float* ptrGlobalResult = result.getPointer();

    if (this->dataset_->getNrows() % 128 != 0 || source_size != this->dataset_->getNrows())
      {
    	throw base::operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
      }

    double time = myOCLKernels->multTransModSPOCL(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims, storageSize);

    // do the rest...
    size_t numWGs = storageSize/OCL_MULT_N_DATAPREFETCH_BLOCKSIZE_SP;
    size_t global = numWGs*OCL_MULT_N_DATAPREFETCH_BLOCKSIZE_SP;

    if (global == 0)
      {
    	global = storageSize;
      }

#pragma omp parallel for
    for (size_t j = global; j < storageSize; j++)
      {
        ptrGlobalResult[j] = 0.0f;

        for (size_t i = 0; i < source_size; i++)
          {
            float curSupport = ptrSource[i];

            for (size_t d = 0; d < dims; d++)
              {
                if (ptrLevel[(j*dims)+d] == 2.0f){
                  curSupport *= 1.0f;
                }
                else if (ptrIndex[(j*dims)+d] == 1.0f){
                  curSupport *= std::max<float>(2.0f - ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d])), 0.0f);
                }
                else if (ptrIndex[(j*dims)+d] == (ptrLevel[(j*dims)+d])-1.0f)
                  {
                    curSupport *= std::max<float>(((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d])) - ptrIndex[(j*dims)+d] + 1.0f, 0.0f);
                  }
                else
                  {
                    curSupport *= std::max<float>(1.0f - fabs( ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d])) - ptrIndex[(j*dims)+d] ), 0.0f);
                  }
              }

            ptrGlobalResult[j] += curSupport;
          }
      }

    return time;
  }

  double OperationMultipleEvalIterativeSPOCLModLinear::multVectorized(base::DataVectorSP& alpha, base::DataVectorSP& result)
  {
    size_t result_size = result.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    result.setAll(0.0f);

    float* ptrAlpha = alpha.getPointer();
    float* ptrData = this->dataset_->getPointer();
    float* ptrResult = result.getPointer();
    float* ptrLevel = this->level_->getPointer();
    float* ptrIndex = this->index_->getPointer();

    if (this->dataset_->getNrows() % 128 != 0 || result_size != this->dataset_->getNrows())
      {
    	throw base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
      }

    double time = myOCLKernels->multModSPOCL(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims, result_size);

    return time;
  }

}

}
