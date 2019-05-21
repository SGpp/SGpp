// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitLinear.hpp>

#include <sgpp/globaldef.hpp>

#include <string.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace sgpp {
namespace pde {

OperationMatrixLTwoDotExplicitLinear::OperationMatrixLTwoDotExplicitLinear() : ownsMatrix_(false) {
  m_ = nullptr;
}

OperationMatrixLTwoDotExplicitLinear::OperationMatrixLTwoDotExplicitLinear(
    sgpp::base::DataMatrix* m, sgpp::base::Grid* grid)
    : ownsMatrix_(false) {
  m_ = m;
  buildMatrix(grid);
}

OperationMatrixLTwoDotExplicitLinear::OperationMatrixLTwoDotExplicitLinear(sgpp::base::Grid* grid)
    : ownsMatrix_(true) {
  m_ = new sgpp::base::DataMatrix(grid->getSize(), grid->getSize());
  buildMatrix(grid);
}

void OperationMatrixLTwoDotExplicitLinear::buildMatrix(sgpp::base::Grid* grid) {
  this->buildMatrixWithBounds(this->m_, grid);
}

OperationMatrixLTwoDotExplicitLinear::~OperationMatrixLTwoDotExplicitLinear() {
  if (ownsMatrix_) delete m_;
}

void OperationMatrixLTwoDotExplicitLinear::mult(sgpp::base::DataVector& alpha,
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
