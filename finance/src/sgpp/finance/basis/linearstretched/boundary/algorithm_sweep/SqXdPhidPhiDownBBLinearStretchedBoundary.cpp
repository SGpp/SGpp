// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/SqXdPhidPhiDownBBLinearStretchedBoundary.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

SqXdPhidPhiDownBBLinearStretchedBoundary::SqXdPhidPhiDownBBLinearStretchedBoundary(
    sgpp::base::GridStorage* storage)
    : SqXdPhidPhiDownBBLinearStretched(storage) {}

SqXdPhidPhiDownBBLinearStretchedBoundary::~SqXdPhidPhiDownBBLinearStretchedBoundary() {}

void SqXdPhidPhiDownBBLinearStretchedBoundary::operator()(sgpp::base::DataVector& source,
                                                          sgpp::base::DataVector& result,
                                                          grid_iterator& index, size_t dim) {
  double q = this->stretching->getIntervalWidth(dim);
  double t = this->stretching->getIntervalOffset(dim);

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

  double bbFactor = ((q * q) + (3.0 * q * t) + (3.0 * t * t)) / (q);

  // check boundary conditions
  if (this->stretching->hasDirichletBoundaryLeft(dim)) {
    result[seq_left] = 0.0;  // left_boundary;
  } else {
    result[seq_left] = (1.0 / 3.0) * left_boundary * bbFactor;
  }

  if (this->stretching->hasDirichletBoundaryRight(dim)) {
    result[seq_right] = 0.0;  // right_boundary;
  } else {
    result[seq_right] = (1.0 / 3.0) * right_boundary * bbFactor;
    // down
    //////////////////////////////////////
    result[seq_right] -= (1.0 / 3.0) * left_boundary * bbFactor;
  }

  // move to root
  if (!index.hint()) {
    index.resetToLevelOne(dim);

    if (!this->storage->isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, left_boundary, right_boundary);
    }

    index.resetToLeftLevelZero(dim);
  }
}

}  // namespace finance
}  // namespace sgpp
