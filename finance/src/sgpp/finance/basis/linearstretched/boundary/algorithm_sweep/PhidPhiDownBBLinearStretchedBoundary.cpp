// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/PhidPhiDownBBLinearStretchedBoundary.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace finance {

PhidPhiDownBBLinearStretchedBoundary::PhidPhiDownBBLinearStretchedBoundary(
    SGPP::base::GridStorage* storage)
    : PhidPhiDownBBLinearStretched(storage) {}

PhidPhiDownBBLinearStretchedBoundary::~PhidPhiDownBBLinearStretchedBoundary() {}

void PhidPhiDownBBLinearStretchedBoundary::operator()(SGPP::base::DataVector& source,
                                                      SGPP::base::DataVector& result,
                                                      grid_iterator& index, size_t dim) {
  // get boundary values
  float_t left_boundary;
  float_t right_boundary;
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

  // the following code is independent from a bounding box

  // check boundary conditions
  if (this->stretching->hasDirichletBoundaryLeft(dim)) {
    result[seq_left] = 0.0;  // left_boundary;
  } else {
    result[seq_left] = left_boundary * (-0.5);
  }

  if (this->stretching->hasDirichletBoundaryRight(dim)) {
    result[seq_right] = 0.0;  // right_boundary;
  } else {
    result[seq_right] = right_boundary * (0.5);
    // down
    //////////////////////////////////////
    result[seq_right] += left_boundary * (0.5);
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

}  // namespace finance
}  // namespace SGPP
