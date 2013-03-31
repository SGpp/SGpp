/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeOCLModLinear.hpp"
#include "base/exception/operation_exception.hpp"

namespace sg
{

namespace parallel
{

  OperationMultipleEvalIterativeOCLModLinear::OperationMultipleEvalIterativeOCLModLinear(
		  base::GridStorage* storage, base::DataMatrix* dataset) :
	  sg::parallel::OperationMultipleEvalVectorized(storage, dataset)
  {
    this->level_ = new base::DataMatrix(storage->size(), storage->dim());
    this->index_ = new base::DataMatrix(storage->size(), storage->dim());

    storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myOCLKernels = new OCLKernels();
  }

  OperationMultipleEvalIterativeOCLModLinear::~OperationMultipleEvalIterativeOCLModLinear()
  {
    delete myOCLKernels;
  }

  void OperationMultipleEvalIterativeOCLModLinear::rebuildLevelAndIndex()
  {
    delete this->level_;
    delete this->index_;

	this->level_ = new base::DataMatrix(storage_->size(), storage_->dim());
	this->index_ = new base::DataMatrix(storage_->size(), storage_->dim());

	storage_->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

    myOCLKernels->resetKernels();
  }

  double OperationMultipleEvalIterativeOCLModLinear::multTransposeVectorized(base::DataVector& source, base::DataVector& result)
  {
    size_t source_size = source.getSize();
	size_t dims = storage_->dim();
	size_t storageSize = storage_->size();

    result.setAll(0.0);

    double* ptrSource = source.getPointer();
    double* ptrData = this->dataset_->getPointer();
    double* ptrLevel = this->level_->getPointer();
    double* ptrIndex = this->index_->getPointer();
    double* ptrGlobalResult = result.getPointer();

    if (this->dataset_->getNrows() % OCL_SGPP_LOCAL_WORKGROUP_SIZE != 0 || source_size != this->dataset_->getNrows())
      {
    	throw base::operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
      }

	size_t numWGs = storageSize/OCL_SGPP_LOCAL_WORKGROUP_SIZE;
    size_t global = numWGs*OCL_SGPP_LOCAL_WORKGROUP_SIZE;

    double time = 0.0;
    if (global > 0)
    	time = myOCLKernels->multTransModOCL(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims, global);

#pragma omp parallel for
    for (size_t j = global; j < storageSize; j++)
      {
        ptrGlobalResult[j] = 0.0f;

        for (size_t i = 0; i < source_size; i++)
          {
            double curSupport = ptrSource[i];

            for (size_t d = 0; d < dims; d++)
              {
                if (ptrLevel[(j*dims)+d] == 2.0)
                  {
                    curSupport *= 1.0;
                  }
                else if (ptrIndex[(j*dims)+d] == 1.0)
                  {
                    curSupport *= std::max<double>(2.0 - ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d])), 0.0);
                  }
                else if (ptrIndex[(j*dims)+d] == (ptrLevel[(j*dims)+d]-1.0))
                  {
                    curSupport *= std::max<double>(((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d])) - ptrIndex[(j*dims)+d] + 1.0, 0.0);
                  }
                else
                  {
                    curSupport *= std::max<double>(1.0 - fabs( ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d])) - ptrIndex[(j*dims)+d] ), 0.0);
                  }
              }

            ptrGlobalResult[j] += curSupport;
          }
      }

    return time;
  }

  double OperationMultipleEvalIterativeOCLModLinear::multVectorized(base::DataVector& alpha, base::DataVector& result)
  {
    size_t result_size = result.getSize();
	size_t dims = storage_->dim();
	size_t storageSize = storage_->size();

    result.setAll(0.0);

    double* ptrAlpha = alpha.getPointer();
    double* ptrData = this->dataset_->getPointer();
    double* ptrResult = result.getPointer();
    double* ptrLevel = this->level_->getPointer();
    double* ptrIndex = this->index_->getPointer();

    if (this->dataset_->getNrows() % OCL_SGPP_LOCAL_WORKGROUP_SIZE != 0 || result_size != this->dataset_->getNrows())
      {
    	throw base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
      }

    double time = myOCLKernels->multModOCL(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims, result_size);

    return time;
  }

}

}
