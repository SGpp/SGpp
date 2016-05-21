// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/XPhidPhiDownBBLinearStretched.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

XPhidPhiDownBBLinearStretched::XPhidPhiDownBBLinearStretched(sgpp::base::GridStorage* storage)
    : storage(storage), stretching(storage->getStretching()) {}

XPhidPhiDownBBLinearStretched::~XPhidPhiDownBBLinearStretched() {}

void XPhidPhiDownBBLinearStretched::operator()(sgpp::base::DataVector& source,
                                               sgpp::base::DataVector& result, grid_iterator& index,
                                               size_t dim) {
  rec(source, result, index, dim, 0.0, 0.0);
}

void XPhidPhiDownBBLinearStretched::rec(sgpp::base::DataVector& source,
                                        sgpp::base::DataVector& result, grid_iterator& index,
                                        size_t dim, double fl, double fr) {
  size_t seq = index.seq();

  double alpha_value = source[seq];

  sgpp::base::GridStorage::index_type::level_type l;
  sgpp::base::GridStorage::index_type::index_type i;

  index.get(dim, l, i);
  // get the positions of the current index as well as its left and right neighbors
  double posl = 0, posr = 0, posc = 0;
  this->stretching->getAdjacentPositions(static_cast<int>(l), static_cast<int>(i), dim, posc, posl,
                                         posr);
  double baseLength = posr - posl;
  double leftLength = posc - posl;

  result[seq] = fl * (1.0 / 6.0) * (2 * posc + 2 * posl - posr) -
                fr * (1.0 / 6.0) * (2 * posc - posl + 2 * posr) -
                1.0 / 6.0 * (baseLength)*alpha_value;  // diagonal entry

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

}  // namespace finance
}  // namespace sgpp
