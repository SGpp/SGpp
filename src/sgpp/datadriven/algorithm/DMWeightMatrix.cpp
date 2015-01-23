/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Zhongwen Song (songz@in.tum.de)
// @author Benjamin (pehersto@in.tum.de)

#include "datadriven/algorithm/DMWeightMatrix.hpp"
#include "base/exception/operation_exception.hpp"
#include "base/operation/BaseOpFactory.hpp"

namespace sg {
  namespace datadriven {

    DMWeightMatrix::DMWeightMatrix(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, sg::base::OperationMatrix& C, double lambda, sg::base::DataVector& w) {
      // create the operations needed in ApplyMatrix
      this->C = &C;
      this->lamb = lambda;
      this->data = &trainData;
      //this->B = SparseGrid.createOperationMultipleEval(this->data);
      this->B = sg::op_factory::createOperationMultipleEval(SparseGrid, *(this->data));
      this->weight = &w;
    }

    DMWeightMatrix::~DMWeightMatrix() {
      delete this->B;
    }


    void DMWeightMatrix::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp((*data).getNrows());
      //size_t M = (*data).getNrows();
      //// Operation B
      this->B->mult(alpha, temp);
      temp.componentwise_mult(*weight);

      this->B->multTranspose(temp, result);

      sg::base::DataVector temptwo(alpha.getSize());
      this->C->mult(alpha, temptwo);
      result.axpy(this->lamb, temptwo);
    }

    void DMWeightMatrix::generateb(sg::base::DataVector& classes, sg::base::DataVector& b) {
      sg::base::DataVector myClassesWithWeights(classes);
      myClassesWithWeights.componentwise_mult(*weight);
      this->B->multTranspose(myClassesWithWeights, b);
    }

  }
}
