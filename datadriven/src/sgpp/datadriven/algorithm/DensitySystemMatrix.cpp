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

DensitySystemMatrix::DensitySystemMatrix(sgpp::base::OperationMatrix* A,
                                         sgpp::base::OperationMultipleEval* B,
                                         sgpp::base::OperationMatrix* C, double lambda,
                                         size_t numSamples)
    : A(A), B(B), C(C), lambda(lambda), numSamples(numSamples) {}

DensitySystemMatrix::DensitySystemMatrix(sgpp::base::Grid& grid, sgpp::base::DataMatrix& trainData,
                                         sgpp::base::OperationMatrix* pC, double lambda)
    : lambda(lambda), numSamples(trainData.getNrows()) {
  A.reset(op_factory::createOperationLTwoDotProduct(grid));
  B.reset(op_factory::createOperationMultipleEval(grid, trainData));
  C.reset(pC);
}

DensitySystemMatrix::~DensitySystemMatrix() {}

void DensitySystemMatrix::mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  result.setAll(0.0);

  // A * alpha
  A->mult(alpha, result);

  // C * alpha
  base::DataVector tmp(result.getSize());
  C->mult(alpha, tmp);

  // A * alpha + lambda * C * alpha
  result.axpy(lambda, tmp);
}

// Matrix-Multiplikation verwenden
void DensitySystemMatrix::generateb(sgpp::base::DataVector& rhs) {
  computeUnweightedRhs(rhs);

  // 1 / M * Bt * 1
  rhs.mult(1. / static_cast<double>(numSamples));
}

void DensitySystemMatrix::computeUnweightedRhs(sgpp::base::DataVector& b) {
  sgpp::base::DataVector y(numSamples);
  y.setAll(1.0);
  // Bt * 1
  B->multTranspose(y, b);
}

}  // namespace datadriven
}  // namespace sgpp
