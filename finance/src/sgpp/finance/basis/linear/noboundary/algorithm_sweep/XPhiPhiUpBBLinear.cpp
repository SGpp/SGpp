// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

XPhiPhiUpBBLinear::XPhiPhiUpBBLinear(sgpp::base::GridStorage* storage)
    : storage(storage), boundingBox(storage->getBoundingBox()) {}

XPhiPhiUpBBLinear::~XPhiPhiUpBBLinear() {}

void XPhiPhiUpBBLinear::operator()(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                   grid_iterator& index, size_t dim) {
  double q = boundingBox->getIntervalWidth(dim);
  double t = boundingBox->getIntervalOffset(dim);

  bool useBB = false;

  if (q != 1.0 || t != 0.0) {
    useBB = true;
  }

  // get boundary values
  double fl = 0.0;
  double fr = 0.0;

  if (useBB) {
    recBB(source, result, index, dim, fl, fr, q, t);
  } else {
    rec(source, result, index, dim, fl, fr);
  }
}

void XPhiPhiUpBBLinear::rec(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                            grid_iterator& index, size_t dim, double& fl, double& fr) {
  size_t seq = index.seq();

  fl = fr = 0.0;
  double fml = 0.0;
  double fmr = 0.0;

  sgpp::base::level_t current_level;
  sgpp::base::index_t current_index;

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

  index.get(dim, current_level, current_index);

  double fm = fml + fmr;

  double alpha_value = source[seq];

  double hsquare = (1.0 / static_cast<double>(1 << (2 * static_cast<int>(current_level))));

  // transposed operations:
  result[seq] = fm;

  fl = ((fm / 2.0) +
        (alpha_value * (hsquare / 2.0 * static_cast<double>(current_index) - hsquare / 12.0))) +
       fl;
  fr = ((fm / 2.0) +
        (alpha_value * (hsquare / 2.0 * static_cast<double>(current_index) + hsquare / 12.0))) +
       fr;
}

void XPhiPhiUpBBLinear::recBB(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                              grid_iterator& index, size_t dim, double& fl, double& fr, double q,
                              double t) {
  size_t seq = index.seq();

  fl = fr = 0.0;
  double fml = 0.0;
  double fmr = 0.0;

  sgpp::base::level_t current_level;
  sgpp::base::index_t current_index;

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      recBB(source, result, index, dim, fl, fml, q, t);
    }

    index.stepRight(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      recBB(source, result, index, dim, fmr, fr, q, t);
    }

    index.up(dim);
  }

  index.get(dim, current_level, current_index);

  double fm = fml + fmr;

  double alpha_value = source[seq];
  double h = (1.0 / static_cast<double>(1 << current_level));

  // transposed operations:
  result[seq] = fm;

  fl = ((fm / 2.0) + (alpha_value * (h * h / 2.0 * static_cast<double>(current_index) * q * q +
                                     h * t * q / 4.0 - h * h / 12.0 * q * q))) +
       fl;
  fr = ((fm / 2.0) + (alpha_value * (h * h / 2.0 * static_cast<double>(current_index) * q * q +
                                     h * t * q / 4.0 + h * h / 12.0 * q * q))) +
       fr;
}

}  // namespace finance
}  // namespace sgpp
