// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/SqXdPhidPhiDownBBLinearStretched.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

SqXdPhidPhiDownBBLinearStretched::SqXdPhidPhiDownBBLinearStretched(sgpp::base::GridStorage* storage)
    : storage(storage), stretching(storage->getStretching()) {}

SqXdPhidPhiDownBBLinearStretched::~SqXdPhidPhiDownBBLinearStretched() {}

void SqXdPhidPhiDownBBLinearStretched::operator()(sgpp::base::DataVector& source,
                                                  sgpp::base::DataVector& result,
                                                  grid_iterator& index, size_t dim) {
  rec(source, result, index, dim, 0.0, 0.0);
}

void SqXdPhidPhiDownBBLinearStretched::rec(sgpp::base::DataVector& source,
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
  double baseLength = posr - posl;
  double leftLength = posc - posl;
  double rightLength = posr - posc;

  double c = 1.0 / 3.0 * (posc + posr + posl);

  result[seq] = 1.0 / 3.0 * alpha_value * baseLength *
                    (posc * (2 * posc + posr + posl) - posr * posl) / (leftLength * rightLength) +
                c * (fl - fr);

  // dehierarchisation
  double fm = (fr - fl) * (leftLength) / (baseLength) + fl + alpha_value;

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
