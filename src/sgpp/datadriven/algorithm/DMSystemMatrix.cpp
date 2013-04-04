/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "datadriven/algorithm/DMSystemMatrix.hpp"
#include "base/exception/operation_exception.hpp"
#include "base/operation/BaseOpFactory.hpp"

namespace sg {
  namespace datadriven {

    DMSystemMatrix::DMSystemMatrix(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, sg::base::OperationMatrix& C, double lambda)
      : DMSystemMatrixBase(trainData, lambda) {
      // create the operations needed in ApplyMatrix
      this->C = &C;
      this->B = sg::op_factory::createOperationMultipleEval(SparseGrid, this->dataset_);
    }

    DMSystemMatrix::~DMSystemMatrix() {
      delete this->B;
    }

    void DMSystemMatrix::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp((*this->dataset_).getNrows());
      size_t M = (*this->dataset_).getNrows();

      // Operation B
      this->B->mult(alpha, temp);
      this->B->multTranspose(temp, result);

      sg::base::DataVector temptwo(alpha.getSize());
      this->C->mult(alpha, temptwo);
      result.axpy(static_cast<double>(M)*this->lambda_, temptwo);
    }

    void DMSystemMatrix::generateb(sg::base::DataVector& classes, sg::base::DataVector& b) {
      this->B->multTranspose(classes, b);
    }

  }

}
