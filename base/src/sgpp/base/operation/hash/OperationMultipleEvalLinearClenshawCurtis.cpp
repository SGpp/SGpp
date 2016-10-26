// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>

#include <sgpp/globaldef.hpp>
#include "OperationMultipleEvalLinearClenshawCurtis.hpp"

namespace sgpp {
namespace base {

void OperationMultipleEvalLinearClenshawCurtis::mult(DataVector& alpha, DataVector& result) {
  AlgorithmDGEMV<SLinearClenshawCurtisBase> op;

  op.mult(storage, base, alpha, this->dataset, result);
}

void OperationMultipleEvalLinearClenshawCurtis::multTranspose(DataVector& source,
                                                              DataVector& result) {
  AlgorithmDGEMV<SLinearClenshawCurtisBase> op;

  op.mult_transposed(storage, base, source, this->dataset, result);
}

double OperationMultipleEvalLinearClenshawCurtis::getDuration() { return 0.0; }

}  // namespace base
}  // namespace sgpp
