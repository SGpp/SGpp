// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/modlinear/algorithm_sweep/PhiPhiDownModLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

PhiPhiDownModLinear::PhiPhiDownModLinear(sgpp::base::GridStorage* storage) : storage(storage) {}

PhiPhiDownModLinear::~PhiPhiDownModLinear() {}

void PhiPhiDownModLinear::operator()(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                     grid_iterator& index, size_t dim) {
  rec(source, result, index, dim, 0.0, 0.0);
}

void PhiPhiDownModLinear::rec(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                              grid_iterator& index, size_t dim, double fl, double fr) {
  size_t seq = index.seq();

  double alpha_value = source[seq];

  sgpp::base::level_t l;
  sgpp::base::index_t i;

  index.get(dim, l, i);

  double h = 1 / pow(2.0, static_cast<int>(l));
  double fm;

  // level 1, constant function
  if (l == 1) {
    // integration
    result[seq] = 0.0 + alpha_value;

    // dehierarchisation
    fm = (fl + fr) / 2.0 + alpha_value;

    // boundary value
    fl += alpha_value;
    fr += alpha_value;
  } else if (i == 1) {  // left boundary
    // integration
    result[seq] = 2.0 / 3.0 * h * (2.0 * fl + fr) + 8.0 / 3.0 * h * alpha_value;

    // dehierarchisation
    fm = (fl + fr) / 2.0 + alpha_value;

    // boundary value
    fl += 2.0 * alpha_value;
  } else if (static_cast<int>(i) == static_cast<int>((1 << l) - 1)) {  // right boundary
    // integration
    result[seq] = 2.0 / 3.0 * h * (fl + 2.0 * fr) + 8.0 / 3.0 * h * alpha_value;

    // dehierarchisation
    fm = (fl + fr) / 2.0 + alpha_value;

    // boundary value
    fr += 2.0 * alpha_value;
  } else {  // inner functions
    // integration
    result[seq] = h * (fl + fr) / 2.0 + 2.0 / 3.0 * h * alpha_value;

    // dehierarchisation
    fm = (fl + fr) / 2.0 + alpha_value;

    // boundary value
  }

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
}  // namespace pde
}  // namespace sgpp
