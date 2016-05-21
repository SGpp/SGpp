// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/modlinear/algorithm_sweep/dPhidPhiUpModLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

dPhidPhiUpModLinear::dPhidPhiUpModLinear(sgpp::base::GridStorage* storage) : storage(storage) {}

dPhidPhiUpModLinear::~dPhidPhiUpModLinear() {}

void dPhidPhiUpModLinear::operator()(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                     grid_iterator& index, size_t dim) {
  double f = 0.0;
  rec(source, result, index, dim, f);
}

void dPhidPhiUpModLinear::rec(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                              grid_iterator& index, size_t dim, double& f) {
  size_t seq = index.seq();

  sgpp::base::GridStorage::index_type::level_type l;
  sgpp::base::GridStorage::index_type::index_type i;

  index.get(dim, l, i);

  double alpha_value = source[seq];
  double ht = pow(2.0, static_cast<int>(l));

  if (l == 1) {
    f = 0.0;

    if (!index.hint()) {
      index.leftChild(dim);

      if (!storage->isValidSequenceNumber(index.seq())) {
        rec(source, result, index, dim, f);
      }

      f = 0.0;
      index.stepRight(dim);

      if (!storage->isValidSequenceNumber(index.seq())) {
        rec(source, result, index, dim, f);
      }

      index.up(dim);
    }

    result[seq] = 0.0;
  } else if (i == 1) {  // left boundary
    f = 0.0;

    if (!index.hint()) {
      index.leftChild(dim);

      if (!storage->isValidSequenceNumber(index.seq())) {
        rec(source, result, index, dim, f);
      }

      index.up(dim);
    }

    result[seq] = ht * f;

    f += 2.0 * alpha_value;
  } else if (static_cast<int>(i) == static_cast<int>((1 << l) - 1)) {  // right boundary
    f = 0.0;

    if (!index.hint()) {
      index.rightChild(dim);

      if (!storage->isValidSequenceNumber(index.seq())) {
        rec(source, result, index, dim, f);
      }

      index.up(dim);
    }

    result[seq] = ht * f;

    f += 2.0 * alpha_value;
  }
}
}  // namespace pde
}  // namespace sgpp
