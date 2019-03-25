/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * OperationMultipleEvalModLinearDistributed.cpp
 *
 * Created on: Mar 24, 2019
 *     Author: Jan Schopohl
 */

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalScalapack/OperationMultipleEvalModLinearDistributed.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/datadriven/algorithm/AlgorithmMultipleEvaluationDistributed.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataVector;
using sgpp::base::LinearModifiedBasis;
using sgpp::base::SLinearModifiedBase;

void OperationMultipleEvalModLinearDistributed::multDistributed(DataVector& alpha,
                                                                DataVectorDistributed& result) {
  // TODO(jan) AlgorithmDGEMV needed?
  AlgorithmMultipleEvaluationDistributed<SLinearModifiedBase> op;
  LinearModifiedBasis<unsigned int, unsigned int> basis;

  op.mult(this->storage, basis, alpha, this->dataset, result);
}

void OperationMultipleEvalModLinearDistributed::multTransposeDistributed(
    DataVector& alpha, DataVectorDistributed& result) {
  // TODO(jan) AlgorithmDGEMV needed?
  AlgorithmMultipleEvaluationDistributed<SLinearModifiedBase> op;
  LinearModifiedBasis<unsigned int, unsigned int> basis;

  op.mult_transpose(this->storage, basis, alpha, this->dataset, result);
}

double OperationMultipleEvalModLinearDistributed::getDuration() { return 0.0; }

}  // namespace datadriven
}  // namespace sgpp