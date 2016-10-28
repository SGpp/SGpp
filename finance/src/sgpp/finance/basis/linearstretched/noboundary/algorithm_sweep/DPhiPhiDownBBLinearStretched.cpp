// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/DPhiPhiDownBBLinearStretched.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

DPhiPhiDownBBLinearStretched::DPhiPhiDownBBLinearStretched(sgpp::base::GridStorage* storage)
    : storage(storage), stretching(storage->getStretching()) {}

DPhiPhiDownBBLinearStretched::~DPhiPhiDownBBLinearStretched() {}

void DPhiPhiDownBBLinearStretched::operator()(sgpp::base::DataVector& source,
                                              sgpp::base::DataVector& result, grid_iterator& index,
                                              size_t dim) {
  rec(source, result, index, dim, 0.0, 0.0);
}

void DPhiPhiDownBBLinearStretched::rec(sgpp::base::DataVector& source,
                                       sgpp::base::DataVector& result, grid_iterator& index,
                                       size_t dim, double fl, double fr) {
  size_t seq = index.seq();

  double alpha_value = source[seq];

  sgpp::base::level_t l;
  sgpp::base::index_t i;

  index.get(dim, l, i);

  double posl = 0, posr = 0, posc = 0;

  this->stretching->getAdjacentPositions(static_cast<int>(l), static_cast<int>(i), dim, posc, posl,
                                         posr);

  // integration
  result[seq] = (0.5 * (fr - fl));  // diagonal entry = 0.0

  // dehierarchisation

  double fm = (fr - fl) * (posc - posl) / (posr - posl) + fl + alpha_value;

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fl, fm);
    }

    index.stepRight(dim);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fm, fr);
    }

    index.up(dim);
  }
}

}  // namespace finance
}  // namespace sgpp
