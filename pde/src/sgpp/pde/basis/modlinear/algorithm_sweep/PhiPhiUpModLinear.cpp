// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/modlinear/algorithm_sweep/PhiPhiUpModLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

PhiPhiUpModLinear::PhiPhiUpModLinear(sgpp::base::GridStorage* storage) : storage(storage) {}

PhiPhiUpModLinear::~PhiPhiUpModLinear() {}

void PhiPhiUpModLinear::operator()(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                   grid_iterator& index, size_t dim) {
  double fl = 0.0;
  double fr = 0.0;
  rec(source, result, index, dim, fl, fr);
}

void PhiPhiUpModLinear::rec(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                            grid_iterator& index, size_t dim, double& fl, double& fr) {
  size_t seq = index.seq();

  double alpha_value = source[seq];

  sgpp::base::level_t l;
  sgpp::base::index_t i;

  index.get(dim, l, i);

  double h = 1 / pow(2.0, static_cast<int>(l));
  double fml = 0.0;
  double fmr = 0.0;

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fl, fml);
    }

    index.stepRight(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fmr, fr);
    }

    index.up(dim);
  }

  double fm = fml + fmr;

  // level 1, constant function
  if (l == 1) {
    result[seq] = fl + fm + fr;

    fl += fm / 2.0 + alpha_value;
    fr += fm / 2.0 + alpha_value;
  } else if (i == 1) {  // left boundary
    result[seq] = 2.0 * fl + fm;

    fl += fm / 2.0 + 4.0 / 3.0 * h * alpha_value;
    fr += fm / 2.0 + 2.0 / 3.0 * h * alpha_value;
  } else if (static_cast<int>(i) == static_cast<int>((1 << l) - 1)) {  // right boundary
    result[seq] = 2.0 * fr + fm;

    fl += fm / 2.0 + 2.0 / 3.0 * h * alpha_value;
    fr += fm / 2.0 + 4.0 / 3.0 * h * alpha_value;
  } else {  // inner functions
    result[seq] = fm;

    fl += fm / 2.0 + h / 2.0 * alpha_value;
    fr += fm / 2.0 + h / 2.0 * alpha_value;
  }
}
}  // namespace pde
}  // namespace sgpp
