// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalScalapack/OperationMultipleEvalKinkLinearDistributed.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearKinkedBasis.hpp>
#include <sgpp/datadriven/algorithm/AlgorithmMultipleEvaluationDistributed.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataVector;
using sgpp::base::LinearKinkedBasis;
using sgpp::base::SLinearKinkedBase;

void OperationMultipleEvalKinkLinearDistributed::multDistributed(DataVector& alpha,
                                                                DataVectorDistributed& result) {
  AlgorithmMultipleEvaluationDistributed<SLinearKinkedBase> op;
  LinearKinkedBasis<unsigned int, unsigned int> basis;

  op.mult(this->storage, basis, alpha, this->dataset, result);
}

void OperationMultipleEvalKinkLinearDistributed::multTransposeDistributed(
    DataVector& alpha, DataVectorDistributed& result) {
  AlgorithmMultipleEvaluationDistributed<SLinearKinkedBase> op;
  LinearKinkedBasis<unsigned int, unsigned int> basis;

  op.mult_transpose(this->storage, basis, alpha, this->dataset, result);
}

double OperationMultipleEvalKinkLinearDistributed::getDuration() { return 0.0; }

}  // namespace datadriven
}  // namespace sgpp
