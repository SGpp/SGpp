// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/SqXdPhidPhiDownBBLinearStretched.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace finance {

SqXdPhidPhiDownBBLinearStretched::SqXdPhidPhiDownBBLinearStretched(SGPP::base::GridStorage* storage)
    : storage(storage), stretching(storage->getStretching()) {}

SqXdPhidPhiDownBBLinearStretched::~SqXdPhidPhiDownBBLinearStretched() {}

void SqXdPhidPhiDownBBLinearStretched::operator()(SGPP::base::DataVector& source,
                                                  SGPP::base::DataVector& result,
                                                  grid_iterator& index, size_t dim) {
  rec(source, result, index, dim, 0.0, 0.0);
}

void SqXdPhidPhiDownBBLinearStretched::rec(SGPP::base::DataVector& source,
                                           SGPP::base::DataVector& result, grid_iterator& index,
                                           size_t dim, float_t fl, float_t fr) {
  size_t seq = index.seq();

  float_t alpha_value = source[seq];

  SGPP::base::GridStorage::index_type::level_type l;
  SGPP::base::GridStorage::index_type::index_type i;

  index.get(dim, l, i);
  float_t posl = 0, posr = 0, posc = 0;

  this->stretching->getAdjacentPositions(static_cast<int>(l), static_cast<int>(i), dim, posc, posl,
                                         posr);
  float_t baseLength = posr - posl;
  float_t leftLength = posc - posl;
  float_t rightLength = posr - posc;

  float_t c = 1.0 / 3.0 * (posc + posr + posl);

  result[seq] = 1.0 / 3.0 * alpha_value * baseLength *
                    (posc * (2 * posc + posr + posl) - posr * posl) / (leftLength * rightLength) +
                c * (fl - fr);

  // dehierarchisation
  float_t fm = (fr - fl) * (leftLength) / (baseLength) + fl + alpha_value;

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->end(index.seq())) {
      rec(source, result, index, dim, fl, fm);
    }

    index.stepRight(dim);

    if (!storage->end(index.seq())) {
      rec(source, result, index, dim, fm, fr);
    }

    index.up(dim);
  }
}

}  // namespace finance
}  // namespace SGPP
