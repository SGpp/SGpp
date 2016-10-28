// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhidPhiUpBBLinearBoundary.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

XPhidPhiUpBBLinearBoundary::XPhidPhiUpBBLinearBoundary(sgpp::base::GridStorage* storage)
    : XPhidPhiUpBBLinear(storage) {}

XPhidPhiUpBBLinearBoundary::~XPhidPhiUpBBLinearBoundary() {}

void XPhidPhiUpBBLinearBoundary::operator()(sgpp::base::DataVector& source,
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

      if (!this->storage->isInvalidSequenceNumber(index.seq())) {
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
      // up
      //////////////////////////////////////
      result[seq_left] = fl;

      result[seq_left] += (source[seq_right] * (((-1.0 / 3.0) * q) - (0.5 * t)));
    }

    if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
      result[seq_right] = 0.0;  // source[seq_right];
    } else {
      result[seq_right] = fr;
    }

    index.resetToLeftLevelZero(dim);
  } else {
    if (!index.hint()) {
      index.resetToLevelOne(dim);

      if (!this->storage->isInvalidSequenceNumber(index.seq())) {
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
      // up
      //////////////////////////////////////
      result[seq_left] = fl;

      result[seq_left] += (source[seq_right] * (-1.0 / 3.0));
    }

    if (this->boundingBox->hasDirichletBoundaryRight(dim)) {
      result[seq_right] = 0.0;  // source[seq_right];
    } else {
      result[seq_right] = fr;
    }

    index.resetToLeftLevelZero(dim);
  }
}

}  // namespace finance
}  // namespace sgpp
