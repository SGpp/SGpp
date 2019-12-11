// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEvalKinkLinear.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearKinkedBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

void OperationMultipleEvalKinkLinear::mult(DataVector& alpha, DataVector& result) {
  AlgorithmDGEMV<SLinearKinkedBase> op;
  LinearKinkedBasis<unsigned int, unsigned int> base;

  op.mult(storage, base, alpha, this->dataset, result);
}

void OperationMultipleEvalKinkLinear::multTranspose(DataVector& source, DataVector& result) {
  AlgorithmDGEMV<SLinearKinkedBase> op;
  LinearKinkedBasis<unsigned int, unsigned int> base;

  op.mult_transposed(storage, base, source, this->dataset, result);
}

double OperationMultipleEvalKinkLinear::getDuration() { return 0.0; }

}  // namespace base
}  // namespace sgpp
