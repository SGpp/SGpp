// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhidPhiUpBBLinearBoundary.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

XdPhidPhiUpBBLinearBoundary::XdPhidPhiUpBBLinearBoundary(sgpp::base::GridStorage* storage)
    : XdPhidPhiUpBBLinear(storage) {}

XdPhidPhiUpBBLinearBoundary::~XdPhidPhiUpBBLinearBoundary() {}

void XdPhidPhiUpBBLinearBoundary::operator()(sgpp::base::DataVector& source,
                                             sgpp::base::DataVector& result, grid_iterator& index,
                                             size_t dim) {
  double q = this->boundingBox->getIntervalWidth(dim);
  double t = this->boundingBox->getIntervalOffset(dim);

  bool useBB = false;

  if (q != 1.0 || t != 0.0) {
    useBB = true;
  }

  // get boundary values
  double fl = 0.0;
  double fr = 0.0;

  if (useBB) {
    if (!index.hint()) {
      index.resetToLevelOne(dim);

      if (!this->storage->isValidSequenceNumber(index.seq())) {
        recBB(source, result, index, dim, fl, fr, q, t);
      }

      index.resetToLeftLevelZero(dim);
    }

    size_t seq_left;
    size_t seq_right;

    // left boundary
    seq_left = index.seq();

    // right boundary
    index.resetToRightLevelZero(dim);
    seq_right = index.seq();

    // check boundary conditions
    if (this->boundingBox->hasDirichletBoundaryLeft(dim)) {
      result[seq_left] = 0.0;  // source[seq_left];
    } else {
      throw base::application_exception(
          "XdPhidPhiUpBBLinearBoundary::operator : Not yet implemented for non-Dirichlet "
          "boundaries.");
    }

    if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
      result[seq_right] = 0.0;  // source[seq_right];
    } else {
      throw base::application_exception(
          "XdPhidPhiUpBBLinearBoundary::operator : Not yet implemented for non-Dirichlet "
          "boundaries.");
    }

    index.resetToLeftLevelZero(dim);
  } else {
    if (!index.hint()) {
      index.resetToLevelOne(dim);

      if (!storage->isValidSequenceNumber(index.seq())) {
        rec(source, result, index, dim, fl, fr);
      }

      index.resetToLeftLevelZero(dim);
    }

    size_t seq_left;
    size_t seq_right;

    // left boundary
    seq_left = index.seq();

    // right boundary
    index.resetToRightLevelZero(dim);
    seq_right = index.seq();

    // check boundary conditions
    if (this->boundingBox->hasDirichletBoundaryLeft(dim)) {
      result[seq_left] = 0.0;  // source[seq_left];
    } else {
      throw base::application_exception(
          "XdPhidPhiUpBBLinearBoundary::operator : Not yet implemented for non-Dirichlet "
          "boundaries.");
    }

    if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
      result[seq_right] = 0.0;  // source[seq_right];
    } else {
      throw base::application_exception(
          "XdPhidPhiUpBBLinearBoundary::operator : Not yet implemented for non-Dirichlet "
          "boundaries.");
    }

    index.resetToLeftLevelZero(dim);
  }
}

}  // namespace finance
}  // namespace sgpp
