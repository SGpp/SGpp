// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitLinearBoundary.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/data_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace pde {

OperationMatrixLTwoDotExplicitLinearBoundary::OperationMatrixLTwoDotExplicitLinearBoundary(
    SGPP::base::DataMatrix* m, SGPP::base::Grid* grid)
    : ownsMatrix_(false) {
  m_ = m;
  buildMatrix(grid);
}

OperationMatrixLTwoDotExplicitLinearBoundary::OperationMatrixLTwoDotExplicitLinearBoundary(
    SGPP::base::Grid* grid)
    : ownsMatrix_(true) {
  m_ = new SGPP::base::DataMatrix(grid->getStorage()->size(), grid->getStorage()->size());
  buildMatrix(grid);
}

OperationMatrixLTwoDotExplicitLinearBoundary::~OperationMatrixLTwoDotExplicitLinearBoundary() {
  if (ownsMatrix_) delete m_;
}

void OperationMatrixLTwoDotExplicitLinearBoundary::buildMatrix(SGPP::base::Grid* grid) {
  // Build matrix (in the moment just by multiplying the OperationMatrix with the unit vectors):
  OperationMatrix* opMatrix = SGPP::op_factory::createOperationLTwoDotProduct(*grid);

  size_t size = grid->getStorage()->size();
  SGPP::base::DataVector unit(size);
  unit.setAll(0.0);
  SGPP::base::DataVector result(size);

  for (size_t i = 0; i < size; i++) {
    // Compute i-th unit vector
    if (i > 0) unit.set(i - 1, 0.0);

    unit.set(i, 1.0);

    // Multiply with operation matrix
    opMatrix->mult(unit, result);
    m_->setColumn(i, result);
  }
}

void OperationMatrixLTwoDotExplicitLinearBoundary::mult(SGPP::base::DataVector& alpha,
                                                        SGPP::base::DataVector& result) {
  size_t nrows = m_->getNrows();
  size_t ncols = m_->getNcols();

  if (alpha.getSize() != ncols || result.getSize() != nrows) {
    throw SGPP::base::data_exception("Dimensions do not match!");
  }

  float_t* data = m_->getPointer();

  // Standard matrix multiplication:
  float_t temp = 0.;
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
}  // namespace SGPP
