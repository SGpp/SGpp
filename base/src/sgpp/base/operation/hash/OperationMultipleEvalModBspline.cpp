// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEvalModBspline.hpp>



#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

void OperationMultipleEvalModBspline::mult(DataVector& alpha,
    DataVector& result) {
  AlgorithmDGEMV<SBsplineModifiedBase> op;

  op.mult(&storage, base, alpha, this->dataset, result);
}

void OperationMultipleEvalModBspline::multTranspose(DataVector& source,
    DataVector& result) {
  AlgorithmDGEMV<SBsplineModifiedBase> op;

  op.mult_transposed(&storage, base, source, this->dataset, result);
}

}  // namespace base
}  // namespace SGPP
