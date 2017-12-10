// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceExplicitLinear.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/DowndPhidPhiBBIterativeLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


#include <sgpp/globaldef.hpp>

#include <string.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace sgpp {
namespace pde {

OperationLaplaceExplicitLinear::OperationLaplaceExplicitLinear(sgpp::base::DataMatrix* m,
                                                               sgpp::base::GridStorage* storage)
  : UpDownOneOpDim(storage), ownsMatrix_(false) {
  m_ = m;
  buildMatrix(storage);
}

OperationLaplaceExplicitLinear::OperationLaplaceExplicitLinear(sgpp::base::GridStorage* storage)
  : UpDownOneOpDim(storage), ownsMatrix_(true) {
  m_ = new sgpp::base::DataMatrix(storage->getSize(), storage->getSize());
  buildMatrix(storage);
}

void OperationLaplaceExplicitLinear::buildMatrix(sgpp::base::GridStorage* storage) {
  size_t ncols = m_->getNcols();
  base::DataVector alpha(ncols);
  base::DataVector beta(ncols);
  // FIXME: inefficient
  for (size_t i = 0; i < ncols; i++) {
    alpha.setAll(0.0);
    alpha.set(i, 1.0);
    UpDownOneOpDim::mult(alpha, beta);
    m_->setColumn(i, beta);
  }
}

OperationLaplaceExplicitLinear::~OperationLaplaceExplicitLinear() {
  if (ownsMatrix_) delete m_;
}

void OperationLaplaceExplicitLinear::mult(sgpp::base::DataVector& alpha,
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

void OperationLaplaceExplicitLinear::specialOP(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result, size_t dim,
                                               size_t gradient_dim) {
  // In direction gradient_dim we only calculate the norm of the gradient
  // The up-part is empty, thus omitted
  if (dim > 0) {
    sgpp::base::DataVector temp(alpha.getSize());
    updown(alpha, temp, dim - 1, gradient_dim);
    downOpDim(temp, result, gradient_dim);
  } else {
    // Terminates dimension recursion
    downOpDim(alpha, result, gradient_dim);
  }
}

void OperationLaplaceExplicitLinear::up(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  PhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<PhiPhiUpBBLinear> s(func, *this->storage);
  s.sweep1D(alpha, result, dim);
}

void OperationLaplaceExplicitLinear::down(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  PhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<PhiPhiDownBBLinear> s(func, *this->storage);
  s.sweep1D(alpha, result, dim);
}

void OperationLaplaceExplicitLinear::downOpDim(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result, size_t dim) {
  DowndPhidPhiBBIterativeLinear myDown(this->storage);
  myDown(alpha, result, dim);
}

void OperationLaplaceExplicitLinear::upOpDim(sgpp::base::DataVector& alpha,
                                             sgpp::base::DataVector& result, size_t dim) {}

}  // namespace pde
}  // namespace sgpp
