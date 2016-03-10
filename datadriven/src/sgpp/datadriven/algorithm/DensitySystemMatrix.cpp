// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinear.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace datadriven {

DensitySystemMatrix::DensitySystemMatrix(sgpp::base::Grid& grid,
    sgpp::base::DataMatrix& trainData, sgpp::base::OperationMatrix& C,
    double lambdaRegression) {
  this->data = &trainData;
  this->lambda = lambdaRegression;

  this->A = sgpp::op_factory::createOperationLTwoDotProduct(grid).release();
  this->B = sgpp::op_factory::createOperationMultipleEval(grid, *(this->data)).release();
  this->C = &C;
}

void DensitySystemMatrix::mult(sgpp::base::DataVector& alpha,
                               sgpp::base::DataVector& result) {
  result.setAll(0.0);

  // A * alpha
  this->A->mult(alpha, result);

  // C * alpha
  base::DataVector tmp(result.getSize());
  this->C->mult(alpha, tmp);

  // A * alpha + lambda * C * alpha
  result.axpy(this->lambda, tmp);
}


// Matrix-Multiplikation verwenden
void DensitySystemMatrix::generateb(sgpp::base::DataVector& rhs) {
  sgpp::base::DataVector y(this->data->getNrows());
  y.setAll(1.0);
  // Bt * 1
  this->B->multTranspose(y, rhs);
  // 1 / 2M * Bt * 1
  rhs.mult(1. / static_cast<double>(this->data->getNrows()));
}

DensitySystemMatrix::~DensitySystemMatrix() {
  delete this->A;
  delete this->B;
}

}  // namespace datadriven
}  // namespace sgpp

