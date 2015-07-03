// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>

#include "OperationMultipleEvalPolyBoundary.hpp"

namespace SGPP {
  namespace base {

    void OperationMultipleEvalPolyBoundary::mult(DataVector& alpha, DataVector& result) {
      AlgorithmDGEMV<SPolyBoundaryBase> op;

      op.mult(storage, base, alpha, this->dataset, result);
    }

    void OperationMultipleEvalPolyBoundary::multTranspose(DataVector& source, DataVector& result) {
      AlgorithmDGEMV<SPolyBoundaryBase> op;

      op.mult_transposed(storage, base, source, this->dataset, result);
    }

  }
}
