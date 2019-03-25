/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * OperationMultipleEvalLinearDistributed.cpp
 *
 * Created on: Mar 23, 2019
 *     Author: Jan Schopohl
 */

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalScalapack/OperationMultipleEvalLinearDistributed.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/datadriven/algorithm/AlgorithmMultipleEvaluationDistributed.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataVector;
using sgpp::base::LinearBasis;
using sgpp::base::SLinearBase;

void OperationMultipleEvalLinearDistributed::multDistributed(DataVector& alpha,
                                                             DataVectorDistributed& result) {
  AlgorithmMultipleEvaluationDistributed<SLinearBase> op;
  LinearBasis<unsigned int, unsigned int> basis;

  op.mult(this->storage, basis, alpha, this->dataset, result);
}

void OperationMultipleEvalLinearDistributed::multTransposeDistributed(
    DataVector& alpha, DataVectorDistributed& result) {
  AlgorithmMultipleEvaluationDistributed<SLinearBase> op;
  LinearBasis<unsigned int, unsigned int> basis;

  op.mult_transpose(this->storage, basis, alpha, this->dataset, result);
}

double OperationMultipleEvalLinearDistributed::getDuration() { return 0.0; }

}  // namespace datadriven
}  // namespace sgpp