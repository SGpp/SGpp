// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/SqrtXPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

SqrtXPhiPhiDownBBLinearBoundary::SqrtXPhiPhiDownBBLinearBoundary(sgpp::base::GridStorage* storage)
    : SqrtXPhiPhiDownBBLinear(storage) {}

SqrtXPhiPhiDownBBLinearBoundary::~SqrtXPhiPhiDownBBLinearBoundary() {}

void SqrtXPhiPhiDownBBLinearBoundary::operator()(sgpp::base::DataVector& source,
                                                 sgpp::base::DataVector& result,
                                                 grid_iterator& index, size_t dim) {
  double q = this->boundingBox->getIntervalWidth(dim);
  double t = this->boundingBox->getIntervalOffset(dim);

  bool useBB = false;

  if (q != 1.0 || t != 0.0) {
    useBB = true;
  }

  // get boundary values
  double left_boundary;
  double right_boundary;
  size_t seq_left;
  size_t seq_right;

  /*
   * Handle Level 0
   */
  // This handles the diagonal only
  //////////////////////////////////////
  // left boundary
  index.resetToLeftLevelZero(dim);
  seq_left = index.seq();
  left_boundary = source[seq_left];

  // right boundary
  index.resetToRightLevelZero(dim);
  seq_right = index.seq();
  right_boundary = source[seq_right];

  if (useBB) {
    // check boundary conditions
    if (this->boundingBox->hasDirichletBoundaryLeft(dim)) {
      result[seq_left] = 0.0;  // left_boundary
    } else {
      throw base::application_exception(
          "SqrtXPhiPhiDownBBLinearBoundary::operator : Not yet implemented for non-Dirichlet "
          "boundaries.");
    }

    if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
      result[seq_right] = 0.0;  // right_boundary;
    } else {
      throw base::application_exception(
          "SqrtXPhiPhiDownBBLinearBoundary::operator : Not yet implemented for non-Dirichlet "
          "boundaries.");
    }

    // move to root
    if (!index.hint()) {
      index.resetToLevelOne(dim);

      if (!this->storage->end(index.seq())) {
        recBB(source, result, index, dim, left_boundary, right_boundary, q, t);
      }

      index.resetToLeftLevelZero(dim);
    }
  } else {
    // check boundary conditions
    if (this->boundingBox->hasDirichletBoundaryLeft(dim)) {
      result[seq_left] = 0.0;  // left_boundary
    } else {
      throw base::application_exception(
          "SqrtXPhiPhiDownBBLinearBoundary::operator : Not yet implemented for non-Dirichlet "
          "boundaries.");
    }

    if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
      result[seq_right] = 0.0;  // right_boundary;
    } else {
      throw base::application_exception(
          "SqrtXPhiPhiDownBBLinearBoundary::operator : Not yet implemented for non-Dirichlet "
          "boundaries.");
    }

    // move to root
    if (!index.hint()) {
      index.resetToLevelOne(dim);

      if (!this->storage->end(index.seq())) {
        rec(source, result, index, dim, left_boundary, right_boundary);
      }

      index.resetToLeftLevelZero(dim);
    }
  }
}

}  // namespace finance
}  // namespace sgpp
