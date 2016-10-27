// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OperationMultipleEvalLinearClenshawCurtisBoundary.hpp"

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>

#include <sgpp/globaldef.hpp>
#include "common/basis/LinearClenshawCurtisBoundaryBasis.hpp"

namespace sgpp {
namespace base {

void OperationMultipleEvalLinearClenshawCurtisBoundary::mult(DataVector& alpha,
                                                             DataVector& result) {
  AlgorithmDGEMV<SLinearClenshawCurtisBoundaryBase> op;

  op.mult(storage, base, alpha, this->dataset, result);
}

void OperationMultipleEvalLinearClenshawCurtisBoundary::multTranspose(DataVector& source,
                                                                      DataVector& result) {
  AlgorithmDGEMV<SLinearClenshawCurtisBoundaryBase> op;

  op.mult_transposed(storage, base, source, this->dataset, result);
}

double OperationMultipleEvalLinearClenshawCurtisBoundary::getDuration() { return 0.0; }

}  // namespace base
}  // namespace sgpp
