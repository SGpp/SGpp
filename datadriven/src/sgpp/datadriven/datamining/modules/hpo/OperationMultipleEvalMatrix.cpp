// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/OperationMultipleEvalMatrix.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

void OperationMultipleEvalMatrix::mult(base::DataVector& alpha, base::DataVector& result) {
  //AlgorithmMultipleEvaluation<SLinearBase> op;
  //LinearBasis<unsigned int, unsigned int> base;

  //op.mult(storage, base, alpha, this->dataset, result);
  this->dataset.mult(alpha, result);
}

void OperationMultipleEvalMatrix::multTranspose(base::DataVector& alpha, base::DataVector& result) {
  //AlgorithmMultipleEvaluation<SLinearBase> op;
  //LinearBasis<unsigned int, unsigned int> base;

  //op.mult_transpose(storage, base, alpha, this->dataset, result);
  sgpp::base::DataMatrix trans{this->dataset};
  trans.transpose();
  trans.mult(alpha, result);
}

double OperationMultipleEvalMatrix::getDuration() { return 0.0; }

}  // namespace datadriven
}  // namespace sgpp
