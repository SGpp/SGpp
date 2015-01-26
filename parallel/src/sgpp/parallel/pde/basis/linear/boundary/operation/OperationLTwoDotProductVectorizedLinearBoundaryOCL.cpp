/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Jacob Jepsen (jepsen@diku.dk)

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/parallel/pde/basis/linear/boundary/operation/OperationLTwoDotProductVectorizedLinearBoundaryOCL.hpp>

#include <cmath>
#include <assert.h>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    OperationLTwoDotProductVectorizedLinearBoundaryOCL::OperationLTwoDotProductVectorizedLinearBoundaryOCL(SGPP::base::GridStorage* storage) : storage(storage) {
      this->OCLPDEKernelsHandle = OCLPDEKernels();
      this->level_ = new SGPP::base::DataMatrix(storage->size(), storage->dim());
      this->level_int_ = new SGPP::base::DataMatrix(storage->size(), storage->dim());
      this->index_ = new SGPP::base::DataMatrix(storage->size(), storage->dim());
      lcl_q = new double[this->storage->dim()];
      storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
      storage->getLevelForIntegral(*(this->level_int_));
    }


    OperationLTwoDotProductVectorizedLinearBoundaryOCL::~OperationLTwoDotProductVectorizedLinearBoundaryOCL() {
      delete this->level_;
      delete this->level_int_;
      delete this->index_;
      delete[] lcl_q;
      this->OCLPDEKernelsHandle.CleanUpGPU();
    }

    void OperationLTwoDotProductVectorizedLinearBoundaryOCL::mult_dirichlet(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      result.setAll(0.0);

      this->OCLPDEKernelsHandle.
      RunOCLKernelLTwoDotBound(alpha, result,
                               this->lcl_q,
                               this->level_->getPointer(),
                               this->index_->getPointer(),
                               this->level_int_->getPointer(),
                               storage->size(),
                               storage->dim(),
                               storage);
    }

    void OperationLTwoDotProductVectorizedLinearBoundaryOCL::mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      result.setAll(0.0);
      bool dirichlet = true;

      // fill q array
      for (size_t d = 0; d < this->storage->dim(); d++) {
        SGPP::base::BoundingBox* boundingBox = this->storage->getBoundingBox();
        lcl_q[d] = boundingBox->getIntervalWidth(d);
        dirichlet = dirichlet && boundingBox->hasDirichletBoundaryLeft(d);
        dirichlet = dirichlet && boundingBox->hasDirichletBoundaryRight(d);
      }

      if (dirichlet) {
        mult_dirichlet(alpha, result);
      } else {
        throw new SGPP::base::operation_exception("OperationLTwoDotProductVectorizedLinearBoundaryOCL::mult : This method is only available on grids with Dirichlet boundaries in all dimensions!");
      }
    }

  }
}

