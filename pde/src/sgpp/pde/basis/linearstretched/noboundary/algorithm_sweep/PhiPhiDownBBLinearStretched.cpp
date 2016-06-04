// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

PhiPhiDownBBLinearStretched::PhiPhiDownBBLinearStretched(sgpp::base::GridStorage* storage)
    : storage(storage), stretching(storage->getStretching()) {}

PhiPhiDownBBLinearStretched::~PhiPhiDownBBLinearStretched() {}

void PhiPhiDownBBLinearStretched::operator()(sgpp::base::DataVector& source,
                                             sgpp::base::DataVector& result, grid_iterator& index,
                                             size_t dim) {
  rec(source, result, index, dim, 0.0, 0.0);
}

void PhiPhiDownBBLinearStretched::rec(sgpp::base::DataVector& source,
                                      sgpp::base::DataVector& result, grid_iterator& index,
                                      size_t dim, double fl, double fr) {
  size_t seq = index.seq();

  double alpha_value = source[seq];

  sgpp::base::level_t current_level;
  sgpp::base::index_t current_index;

  index.get(dim, current_level, current_index);

  double posl = 0, posr = 0, currentPosition = 0;
  this->stretching->getAdjacentPositions(static_cast<int>(current_level),
                                         static_cast<int>(current_index), dim, currentPosition,
                                         posl, posr);
  double baseLength = posr - posl;
  double leftLength = currentPosition - posl;
  double rightLength = posr - currentPosition;

  // integration
  result[seq] = (1.0 / 3.0) * (baseLength)*alpha_value + fl / 6.0 * (baseLength + rightLength) +
                fr / 6.0 * (baseLength + leftLength);

  // dehierarchisation
  double fm = (fr - fl) * (leftLength) / (baseLength) + fl + alpha_value;

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fl, fm);
    }

    index.stepRight(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fm, fr);
    }

    index.up(dim);
  }
}

}  // namespace pde
}  // namespace sgpp
