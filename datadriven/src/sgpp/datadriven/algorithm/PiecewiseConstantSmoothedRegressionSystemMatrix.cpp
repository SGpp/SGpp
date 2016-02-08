// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/PiecewiseConstantSmoothedRegressionSystemMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinear.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

PiecewiseConstantSmoothedRegressionSystemMatrix::PiecewiseConstantSmoothedRegressionSystemMatrix(
  SGPP::datadriven::PiecewiseConstantRegression::Node& piecewiseRegressor,
  SGPP::base::Grid& grid, SGPP::base::OperationMatrix& C,
  float_t lambdaRegression) :
  piecewiseRegressor(piecewiseRegressor), grid(grid) {
  this->lambda = lambdaRegression;

  this->A = SGPP::op_factory::createOperationLTwoDotProduct(grid);
  //      this->B = SGPP::op_factory::createOperationMultipleEval(grid, *(this->data));
  this->C = &C;
}

void PiecewiseConstantSmoothedRegressionSystemMatrix::mult(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  result.setAll(0.0);

  // A * alpha
  this->A->mult(alpha, result);

  // C * alpha
  base::DataVector tmp(result.getSize());
  this->C->mult(alpha, tmp);

  // A * alpha + lambda * C * alpha
  result.axpy(this->lambda, tmp);
}

// Matrix-Multiplikation verwenden
void PiecewiseConstantSmoothedRegressionSystemMatrix::generateb(
  SGPP::base::DataVector& rhs) {

  //store result in rhs!
  SGPP::base::GridStorage* storage = grid.getStorage();
  uint64_t totalIntegratedNodes = 0;
  #pragma omp parallel for

  for (size_t gridIndex = 0; gridIndex < storage->size(); gridIndex++) {
    SGPP::base::GridIndex* gridPoint = storage->get(gridIndex);
    size_t integratedNodes;
    rhs[gridIndex] = piecewiseRegressor.integrate(*gridPoint, integratedNodes);
    #pragma omp atomic
    totalIntegratedNodes += integratedNodes;
  }

  std::cout << "totalIntegratedNodes: " << totalIntegratedNodes << std::endl;
  std::cout << "integrated nodes per grid point: "
            << (static_cast<float_t>(totalIntegratedNodes) / static_cast<float_t>
                (storage->size())) << std::endl;
}

PiecewiseConstantSmoothedRegressionSystemMatrix::~PiecewiseConstantSmoothedRegressionSystemMatrix() {
  delete this->A;
}

}
}
