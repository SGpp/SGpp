// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>

#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPoly.hpp>



#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

void OperationMultipleEvalPoly::mult(DataVector& alpha, DataVector& result) {
  AlgorithmDGEMV<SPolyBase> op;

  op.mult(storage, base, alpha, this->dataset, result);
}

void OperationMultipleEvalPoly::multTranspose(DataVector& source,
    DataVector& result) {
  AlgorithmDGEMV<SPolyBase> op;

  op.mult_transposed(storage, base, source, this->dataset, result);
}

}  // namespace base
}  // namespace sgpp
