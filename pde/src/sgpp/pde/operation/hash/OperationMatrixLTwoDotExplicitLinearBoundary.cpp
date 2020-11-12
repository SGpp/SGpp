// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitLinearBoundary.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/data_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationMatrixLTwoDotExplicitLinearBoundary::OperationMatrixLTwoDotExplicitLinearBoundary(
    sgpp::base::DataMatrix* m, sgpp::base::Grid* grid)
    : ownsMatrix_(false) {
  m_ = m;
  buildMatrix(grid);
}

OperationMatrixLTwoDotExplicitLinearBoundary::OperationMatrixLTwoDotExplicitLinearBoundary(
    sgpp::base::Grid* grid)
    : ownsMatrix_(true) {
  m_ = new sgpp::base::DataMatrix(grid->getSize(), grid->getSize());
  buildMatrix(grid);
}

OperationMatrixLTwoDotExplicitLinearBoundary::~OperationMatrixLTwoDotExplicitLinearBoundary() {
  if (ownsMatrix_) delete m_;
}

void OperationMatrixLTwoDotExplicitLinearBoundary::buildMatrix(sgpp::base::Grid* grid) {
  // Build matrix (in the moment just by multiplying the OperationMatrix with the unit vectors):
  std::unique_ptr<OperationMatrix> opMatrix(
      sgpp::op_factory::createOperationLTwoDotProduct(*grid));

  size_t size = grid->getSize();
  sgpp::base::DataVector unit(size);
  sgpp::base::DataVector result(size);

  for (size_t i = 0; i < size; i++) {
    // Compute i-th unit vector
    if (i > 0) unit.set(i - 1, 0.0);

    unit.set(i, 1.0);

    // Multiply with operation matrix
    opMatrix->mult(unit, result);
    m_->setColumn(i, result);
  }
}

void OperationMatrixLTwoDotExplicitLinearBoundary::mult(sgpp::base::DataVector& alpha,
                                                        sgpp::base::DataVector& result) {
  size_t nrows = m_->getNrows();
  size_t ncols = m_->getNcols();

  if (alpha.getSize() != ncols || result.getSize() != nrows) {
    throw sgpp::base::data_exception("Dimensions do not match!");
  }

  double* data = m_->getPointer();

  // Standard matrix multiplication:
  double temp = 0.;
  size_t acc = 0;

  for (size_t i = 0; i < nrows; i++) {
    for (size_t j = 0; j < ncols; j++) {
      temp += data[j + acc] * alpha[j];
    }

    result[i] = temp;
    temp = 0.;
    acc += ncols;
  }
}

}  // namespace pde
}  // namespace sgpp
