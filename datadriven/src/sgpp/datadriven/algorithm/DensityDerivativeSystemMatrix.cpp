// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DensityDerivativeSystemMatrix.hpp>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinear.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

DensityDerivativeSystemMatrix::DensityDerivativeSystemMatrix(sgpp::base::OperationMatrix* A,
                                                             sgpp::base::OperationMultipleEval* B,
                                                             sgpp::base::OperationMatrix* C,
                                                             double lambda, size_t numSamples)
    : A(A), B(B), C(C), lambda(lambda), numSamples(numSamples) {}

DensityDerivativeSystemMatrix::DensityDerivativeSystemMatrix(sgpp::base::Grid& grid,
                                                             sgpp::base::DataMatrix& trainData,
                                                             sgpp::base::OperationMatrix* pC,
                                                             double lambda, size_t derivDim)
    : lambda(lambda), numSamples(trainData.getNrows()) {
  A.reset(op_factory::createOperationLTwoDotProduct(grid));
  // B is an evaluation of the partial derivative
  B.reset(op_factory::createOperationMultipleEvalPartialDerivativeNaive(grid, trainData, derivDim));
  C.reset(pC);
}

DensityDerivativeSystemMatrix::~DensityDerivativeSystemMatrix() {}

void DensityDerivativeSystemMatrix::mult(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataVector& result) {
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
void DensityDerivativeSystemMatrix::generateb(sgpp::base::DataVector& rhs) {
  // Compute unweighted rhs; contains the sign!!!
  computeUnweightedRhs(rhs);

  // 1 / M * Bt * (-1)
  rhs.mult(1. / static_cast<double>(numSamples));
}

void DensityDerivativeSystemMatrix::computeUnweightedRhs(sgpp::base::DataVector& b) {
  sgpp::base::DataVector y(numSamples);
  // In general: (-1)^|j| / M
  // For first derivative: |j| = 1
  // Bt * (-1)
  y.setAll(-1.0);
  B->multTranspose(y, b);
}

}  // namespace datadriven
}  // namespace sgpp
