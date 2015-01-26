/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Jacob Jepsen (jepsen@diku.dk)

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/parallel/pde/basis/linear/noboundary/operation/OperationLTwoDotProductVectorizedLinearOCL.hpp>

#include <cmath>

namespace sg {
  namespace parallel {

    OperationLTwoDotProductVectorizedLinearOCL::OperationLTwoDotProductVectorizedLinearOCL(sg::base::GridStorage* storage) : storage(storage) {

      this->lcl_q = new double[this->storage->dim()];
      this->OCLPDEKernelsHandle = OCLPDEKernels();

      this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      this->level_int_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());

      storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
      storage->getLevelForIntegral(*(this->level_int_));

    }

    OperationLTwoDotProductVectorizedLinearOCL::~OperationLTwoDotProductVectorizedLinearOCL() {
      delete this->level_;
      delete this->level_int_;
      delete this->index_;
      delete[] lcl_q;
      this->OCLPDEKernelsHandle.CleanUpGPU();
    }

    void OperationLTwoDotProductVectorizedLinearOCL::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);

      for (size_t d = 0; d < this->storage->dim(); d++) {
        sg::base::BoundingBox* boundingBox = this->storage->getBoundingBox();
        this->lcl_q[d] = boundingBox->getIntervalWidth(d);
      }

      this->OCLPDEKernelsHandle.
      RunOCLKernelLTwoDotInner(alpha, result,
                               this->lcl_q,
                               this->level_->getPointer(),
                               this->index_->getPointer(),
                               this->level_int_->getPointer(),
                               storage->size(),
                               storage->dim(),
                               storage);
    }

  }
}
