// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/PhidPhiDownBBLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

PhidPhiDownBBLinear::PhidPhiDownBBLinear(sgpp::base::GridStorage* storage)
    : storage(storage), boundingBox(storage->getBoundingBox()) {}

PhidPhiDownBBLinear::~PhidPhiDownBBLinear() {}

void PhidPhiDownBBLinear::operator()(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                     grid_iterator& index, size_t dim) {
  rec(source, result, index, dim, 0.0, 0.0);
}

void PhidPhiDownBBLinear::rec(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                              grid_iterator& index, size_t dim, double fl, double fr) {
  size_t seq = index.seq();

  double alpha_value = source[seq];

  sgpp::base::level_t l;
  sgpp::base::index_t i;

  index.get(dim, l, i);

  // integration
  result[seq] = (0.5 * (fl - fr));  // diagonal entry = 0.0

  // dehierarchisation
  double fm = ((fl + fr) / 2.0) + alpha_value;

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

}  // namespace finance
}  // namespace sgpp
