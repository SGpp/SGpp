// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/modlinear/algorithm_sweep/dPhidPhiDownModLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

dPhidPhiDownModLinear::dPhidPhiDownModLinear(sgpp::base::GridStorage* storage) : storage(storage) {}

dPhidPhiDownModLinear::~dPhidPhiDownModLinear() {}

void dPhidPhiDownModLinear::operator()(sgpp::base::DataVector& source,
                                       sgpp::base::DataVector& result, grid_iterator& index,
                                       size_t dim) {
  rec(source, result, index, dim, 0.0);
}

void dPhidPhiDownModLinear::rec(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                grid_iterator& index, size_t dim, double f) {
  size_t seq = index.seq();
  sgpp::base::GridStorage::index_type::level_type l;
  sgpp::base::GridStorage::index_type::index_type i;

  index.get(dim, l, i);

  double alpha_value = source[seq];
  double ht = pow(2.0, static_cast<int>(l));
  double f_local = 0.0;

  // level 1, constant function
  if (l == 1) {
    f_local = 0.0;
    result[seq] = 0.0 + 0.0;
  } else if ((i == 1) || (static_cast<int>(i) == static_cast<int>((1 << l) - 1))) {
    // left boundary & right boundary
    f_local = ht * alpha_value;
    result[seq] = 2.0 * f + 2.0 * f_local;
  } else {
    // inner functions
    f_local = ht * alpha_value;
    result[seq] = 0.0 + 2.0 * f_local;
  }

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, f + f_local);
    }

    index.stepRight(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, f + f_local);
    }

    index.up(dim);
  }
}
}  // namespace pde
}  // namespace sgpp
