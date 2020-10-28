// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DensityDifferenceSystemMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinear.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

DensityDifferenceSystemMatrix::DensityDifferenceSystemMatrix(sgpp::base::OperationMatrix* A,
                                                             sgpp::base::OperationMultipleEval* B_p,
                                                             sgpp::base::OperationMultipleEval* B_q,
                                                             sgpp::base::OperationMatrix* C,
                                                             double lambda, size_t numSamples_P,
                                                             size_t numSamples_Q)
    : A(A),
      B_p(B_p),
      B_q(B_q),
      C(C),
      lambda(lambda),
      numSamples_P(numSamples_P),
      numSamples_Q(numSamples_Q) {}

DensityDifferenceSystemMatrix::DensityDifferenceSystemMatrix(sgpp::base::Grid& grid,
                                                             sgpp::base::DataMatrix& trainData_P,
                                                             sgpp::base::DataMatrix& trainData_Q,
                                                             sgpp::base::OperationMatrix* pC,
                                                             double lambda)
    : lambda(lambda), numSamples_P(trainData_P.getNrows()), numSamples_Q(trainData_Q.getNrows()) {
  A.reset(op_factory::createOperationLTwoDotProduct(grid));
  B_p.reset(op_factory::createOperationMultipleEval(grid, trainData_P));
  B_q.reset(op_factory::createOperationMultipleEval(grid, trainData_Q));
  C.reset(pC);
}

DensityDifferenceSystemMatrix::~DensityDifferenceSystemMatrix() {}

void DensityDifferenceSystemMatrix::mult(sgpp::base::DataVector& alpha,
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
void DensityDifferenceSystemMatrix::generateb(sgpp::base::DataVector& rhs) {
  sgpp::base::DataVector rhs2(rhs.getSize());
  computeUnweightedRhs(rhs, rhs2);

  // (1 / M_p) * Bt_p * 1
  rhs.mult(1. / static_cast<double>(numSamples_P));

  // (1 / M_p) * Bt_p * 1 - (1 / M_q) * Bt_q * 1
  rhs.axpy(-1. / static_cast<double>(numSamples_Q), rhs2);
}

void DensityDifferenceSystemMatrix::computeUnweightedRhs(base::DataVector& bp,
                                                         base::DataVector& bq) {
  sgpp::base::DataVector yb(numSamples_P);
  yb.setAll(1.0);
  // Bt_p * 1
  B_p->multTranspose(yb, bp);

  sgpp::base::DataVector yq(numSamples_Q);
  yq.setAll(1.0);
  // Bt_q * 1
  B_q->multTranspose(yq, bq);
}

}  // namespace datadriven
}  // namespace sgpp
