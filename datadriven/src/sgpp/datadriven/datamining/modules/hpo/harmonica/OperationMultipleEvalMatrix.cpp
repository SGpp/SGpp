// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/harmonica/OperationMultipleEvalMatrix.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>

namespace sgpp {
namespace datadriven {

void OperationMultipleEvalMatrix::mult(base::DataVector &alpha, base::DataVector &result) {
  this->dataset.mult(alpha, result);
}

void OperationMultipleEvalMatrix::multTranspose(base::DataVector &alpha, base::DataVector &result) {
  sgpp::base::DataMatrix trans(this->dataset);
  trans.transpose();
  trans.mult(alpha, result);
}

double OperationMultipleEvalMatrix::getDuration() { return 0.0; }
}  // namespace datadriven
}  // namespace sgpp
