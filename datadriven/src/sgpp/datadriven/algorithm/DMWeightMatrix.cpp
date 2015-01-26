/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Zhongwen Song (songz@in.tum.de)
// @author Benjamin (pehersto@in.tum.de)

#include <sgpp/datadriven/algorithm/DMWeightMatrix.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    DMWeightMatrix::DMWeightMatrix(SGPP::base::Grid& SparseGrid, SGPP::base::DataMatrix& trainData, SGPP::base::OperationMatrix& C, double lambda, SGPP::base::DataVector& w) {
      // create the operations needed in ApplyMatrix
      this->C = &C;
      this->lamb = lambda;
      this->data = &trainData;
      //this->B = SparseGrid.createOperationMultipleEval(this->data);
      this->B = SGPP::op_factory::createOperationMultipleEval(SparseGrid, *(this->data));
      this->weight = &w;
    }

    DMWeightMatrix::~DMWeightMatrix() {
      delete this->B;
    }


    void DMWeightMatrix::mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      SGPP::base::DataVector temp((*data).getNrows());
      //size_t M = (*data).getNrows();
      //// Operation B
      this->B->mult(alpha, temp);
      temp.componentwise_mult(*weight);

      this->B->multTranspose(temp, result);

      SGPP::base::DataVector temptwo(alpha.getSize());
      this->C->mult(alpha, temptwo);
      result.axpy(this->lamb, temptwo);
    }

    void DMWeightMatrix::generateb(SGPP::base::DataVector& classes, SGPP::base::DataVector& b) {
      SGPP::base::DataVector myClassesWithWeights(classes);
      myClassesWithWeights.componentwise_mult(*weight);
      this->B->multTranspose(myClassesWithWeights, b);
    }

  }
}
