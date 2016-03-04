// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/LogDensitySystemMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinear.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

LogDensitySystemMatrix::LogDensitySystemMatrix(base::Grid& grid, base::DataVector& alphaRef,
                                               base::DataMatrix& trainData,
                                               base::OperationMatrix& C, double lambdaRegression)
    : grid(grid), alpha(alphaRef), data(trainData), lambda(lambdaRegression), C(C) {
  A = op_factory::createOperationLTwoDotProduct(grid);
  B = op_factory::createOperationMultipleEval(grid, data);
}

LogDensitySystemMatrix::~LogDensitySystemMatrix() {}

void LogDensitySystemMatrix::mult(base::DataVector& alpha, base::DataVector& result) {
  result.setAll(0.0);

  // A * alpha
  A->mult(alpha, result);

  // C * alpha
  base::DataVector tmp(result.getSize());
  C.mult(alpha, tmp);

  // A * alpha + lambda * C * alpha
  result.axpy(lambda, tmp);
}

void LogDensitySystemMatrix::generateb(base::DataVector& rhs) {
  base::DataVector y(data.getNrows());

  // compute y_i = log(p(x_i)) / p(x_i)
  double val = 0.0;
  std::unique_ptr<base::OperationEval> opEval = op_factory::createOperationEval(grid);
  base::DataVector point(data.getNcols());

  for (size_t i = 0; i < y.getSize(); i++) {
    data.getRow(i, point);
    val = opEval->eval(alpha, point);

    if (val < 1e-10) {
      // set negative values manually to -1e10
      y[i] = 0.0;
    } else {
      y[i] = std::log(val) / val;
    }
  }

  // Bt * y
  B->multTranspose(y, rhs);
  // 1 / 2M * Bt * 1
  rhs.mult(1. / static_cast<double>(data.getNrows()));
}
}  // namespace datadriven
}  // namespace sgpp
