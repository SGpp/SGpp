// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinear.hpp>


#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

void OperationMultipleEvalLinear::mult(DataVector& alpha, DataVector& result) {
  AlgorithmMultipleEvaluation<SLinearBase> op;
  LinearBasis<unsigned int, unsigned int> base;

  op.mult(storage, base, alpha, this->dataset, result);
}

void OperationMultipleEvalLinear::multTranspose(DataVector& alpha,
    DataVector& result) {
  AlgorithmMultipleEvaluation<SLinearBase> op;
  LinearBasis<unsigned int, unsigned int> base;

  op.mult_transpose(storage, base, alpha, this->dataset, result);
}

}  // namespace base
}  // namespace sgpp
