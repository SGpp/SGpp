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


namespace SGPP {
namespace datadriven {

DensitySystemMatrix::DensitySystemMatrix(SGPP::base::Grid& grid,
    SGPP::base::DataMatrix& trainData, SGPP::base::OperationMatrix& C,
    float_t lambdaRegression) {
  this->data = &trainData;
  this->lambda = lambdaRegression;

  this->A = SGPP::op_factory::createOperationLTwoDotProduct(grid).release();
  this->B = SGPP::op_factory::createOperationMultipleEval(grid, *(this->data)).release();
  this->C = &C;
}

void DensitySystemMatrix::mult(SGPP::base::DataVector& alpha,
                               SGPP::base::DataVector& result) {
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
void DensitySystemMatrix::generateb(SGPP::base::DataVector& rhs) {
  SGPP::base::DataVector y(this->data->getNrows());
  y.setAll(1.0);
  // Bt * 1
  this->B->multTranspose(y, rhs);
  // 1 / 2M * Bt * 1
  rhs.mult(1. / (float_t)this->data->getNrows());
}

DensitySystemMatrix::~DensitySystemMatrix() {
  delete this->A;
  delete this->B;
}

}  // namespace datadriven
}  // namespace SGPP

