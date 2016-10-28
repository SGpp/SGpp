// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

PhiPhiDownBBLinearBoundary::PhiPhiDownBBLinearBoundary(sgpp::base::GridStorage* storage)
    : PhiPhiDownBBLinear(storage) {}

PhiPhiDownBBLinearBoundary::~PhiPhiDownBBLinearBoundary() {}

void PhiPhiDownBBLinearBoundary::operator()(sgpp::base::DataVector& source,
                                            sgpp::base::DataVector& result, grid_iterator& index,
                                            size_t dim) {
  // std::cout << dim << std::endl;
  // std::cout << index.toString() << std::endl;

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
      result[seq_left] = ((1.0 / 3.0) * left_boundary) * q;
    }

    if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
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

      if (!this->storage->isInvalidSequenceNumber(index.seq())) {
        recBB(source, result, index, dim, left_boundary, right_boundary, q, t);
      }

      index.resetToLeftLevelZero(dim);
    }
  } else {
    // check boundary conditions
    if (this->boundingBox->hasDirichletBoundaryLeft(dim)) {
      result[seq_left] = 0.0;  // left_boundary
    } else {
      result[seq_left] = (1.0 / 3.0) * left_boundary;
    }

    if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
      result[seq_right] = 0.0;  // right_boundary;
    } else {
      result[seq_right] = (1.0 / 3.0) * right_boundary;

      // down
      //////////////////////////////////////
      result[seq_right] += (1.0 / 6.0) * left_boundary;
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
}

}  // namespace pde
}  // namespace sgpp
