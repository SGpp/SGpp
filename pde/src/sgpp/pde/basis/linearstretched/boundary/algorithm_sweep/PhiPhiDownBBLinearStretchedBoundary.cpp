// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

sgpp::pde::PhiPhiDownBBLinearStretchedBoundary::PhiPhiDownBBLinearStretchedBoundary(
    sgpp::base::GridStorage* storage)
    : sgpp::pde::PhiPhiDownBBLinearStretched(storage) {}

sgpp::pde::PhiPhiDownBBLinearStretchedBoundary::~PhiPhiDownBBLinearStretchedBoundary() {}

void sgpp::pde::PhiPhiDownBBLinearStretchedBoundary::operator()(sgpp::base::DataVector& source,
                                                                sgpp::base::DataVector& result,
                                                                grid_iterator& index, size_t dim) {
  double q = this->stretching->getIntervalWidth(dim);
  //  double t = this->stretching->getIntervalOffset(dim);

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

  // check boundary conditions
  if (this->stretching->hasDirichletBoundaryLeft(dim)) {
    result[seq_left] = 0.0;  // left_boundary
  } else {
    result[seq_left] = ((1.0 / 3.0) * left_boundary) * q;
  }

  if (this->stretching->hasDirichletBoundaryRight(dim)) {
    result[seq_right] = 0.0;  // right_boundary;
  } else {
    result[seq_right] = ((1.0 / 3.0) * right_boundary) * q;

    // down
    //////////////////////////////////////
    result[seq_right] += ((1.0 / 6.0) * left_boundary) * q;
  }

  // move to root
  if (!index.hint()) {
    index.resetToLevelOne(dim);

    if (!this->storage->isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, left_boundary, right_boundary);
    }

    index.resetToLeftLevelZero(dim);
  }
}

}  // namespace pde
}  // namespace sgpp
