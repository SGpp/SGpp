/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Jacob Jepsen (jepsen@diku.dk)

#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/generation/GridGenerator.hpp"
#include "base/exception/operation_exception.hpp"

#include "parallel/pde/basis/linear/boundary/operation/OperationLaplaceVectorizedLinearBoundaryOCL.hpp"

#include <cmath>
#include <assert.h>

namespace sg {
  namespace parallel {

    OperationLaplaceVectorizedLinearBoundaryOCL::OperationLaplaceVectorizedLinearBoundaryOCL(sg::base::GridStorage* storage, sg::base::DataVector& lambda) : storage(storage) {

      this->lambda = new sg::base::DataVector(lambda);
      this->OCLPDEKernelsHandle = OCLPDEKernels();
      this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      this->level_int_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      lcl_q = new double[this->storage->dim()];
      lcl_q_inv = new double[this->storage->dim()];
      storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
      storage->getLevelForIntegral(*(this->level_int_));

    }

    OperationLaplaceVectorizedLinearBoundaryOCL::OperationLaplaceVectorizedLinearBoundaryOCL(sg::base::GridStorage* storage) : storage(storage) {


      this->lambda = new base::DataVector(storage->dim());
      this->lambda->setAll(1.0);
      this->OCLPDEKernelsHandle = OCLPDEKernels();
      this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      this->level_int_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      lcl_q = new double[this->storage->dim()];
      lcl_q_inv = new double[this->storage->dim()];
      storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
      storage->getLevelForIntegral(*(this->level_int_));

    }


    OperationLaplaceVectorizedLinearBoundaryOCL::~OperationLaplaceVectorizedLinearBoundaryOCL() {
      delete this->level_;
      delete this->level_int_;
      delete this->index_;
      delete[] lcl_q;
      delete[] lcl_q_inv;

      this->OCLPDEKernelsHandle.CleanUpGPU();
    }

    void OperationLaplaceVectorizedLinearBoundaryOCL::mult_dirichlet(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);

      this->OCLPDEKernelsHandle.RunOCLKernelLaplaceBound(alpha, result, lcl_q, lcl_q_inv,
          this->level_->getPointer(),
          this->index_->getPointer(),
          this->level_int_->getPointer(),
          lambda->getPointer(),
          storage->size(),
          storage->dim(),
          storage);


    }

    void OperationLaplaceVectorizedLinearBoundaryOCL::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);
      bool dirichlet = true;

      // fill q array
      for (size_t d = 0; d < this->storage->dim(); d++) {
        sg::base::BoundingBox* boundingBox = this->storage->getBoundingBox();
        lcl_q[d] = boundingBox->getIntervalWidth(d);
        lcl_q_inv[d] = 1.0 / boundingBox->getIntervalWidth(d);
        dirichlet = dirichlet && boundingBox->hasDirichletBoundaryLeft(d);
        dirichlet = dirichlet && boundingBox->hasDirichletBoundaryRight(d);
      }

      if (dirichlet) {
        mult_dirichlet(alpha, result);
      } else {
        throw new sg::base::operation_exception("OperationLaplaceVectorizedLinearBoundaryOCL::mult : This method is only available on grids with Dirichlet boundaries in all dimensions!");
      }
    }
  }
}
