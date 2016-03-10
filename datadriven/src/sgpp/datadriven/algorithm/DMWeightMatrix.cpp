// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DMWeightMatrix.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace datadriven {

DMWeightMatrix::DMWeightMatrix(sgpp::base::Grid& SparseGrid,
                               sgpp::base::DataMatrix& trainData, sgpp::base::OperationMatrix& C,
                               double lambda, sgpp::base::DataVector& w) {
  // create the operations needed in ApplyMatrix
  this->C = &C;
  this->lamb = lambda;
  this->data = &trainData;
  // this->B = SparseGrid.createOperationMultipleEval(this->data);
  this->B = sgpp::op_factory::createOperationMultipleEval(SparseGrid,
            *(this->data)).release();
  this->weight = &w;
}

DMWeightMatrix::~DMWeightMatrix() {
  delete this->B;
}


void DMWeightMatrix::mult(sgpp::base::DataVector& alpha,
                          sgpp::base::DataVector& result) {
  sgpp::base::DataVector temp((*data).getNrows());
  // size_t M = (*data).getNrows();
  //// Operation B
  this->B->mult(alpha, temp);
  temp.componentwise_mult(*weight);

  this->B->multTranspose(temp, result);

  sgpp::base::DataVector temptwo(alpha.getSize());
  this->C->mult(alpha, temptwo);
  result.axpy(this->lamb, temptwo);
}

void DMWeightMatrix::generateb(sgpp::base::DataVector& classes,
                               sgpp::base::DataVector& b) {
  sgpp::base::DataVector myClassesWithWeights(classes);
  myClassesWithWeights.componentwise_mult(*weight);
  this->B->multTranspose(myClassesWithWeights, b);
}

}  // namespace datadriven
}  // namespace sgpp
