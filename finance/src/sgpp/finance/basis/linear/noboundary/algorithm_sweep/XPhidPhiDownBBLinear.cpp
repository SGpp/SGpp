// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhidPhiDownBBLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

XPhidPhiDownBBLinear::XPhidPhiDownBBLinear(sgpp::base::GridStorage* storage)
    : storage(storage), boundingBox(storage->getBoundingBox()) {}

XPhidPhiDownBBLinear::~XPhidPhiDownBBLinear() {}

void XPhidPhiDownBBLinear::operator()(sgpp::base::DataVector& source,
                                      sgpp::base::DataVector& result, grid_iterator& index,
                                      size_t dim) {
  double q = boundingBox->getIntervalWidth(dim);
  double t = boundingBox->getIntervalOffset(dim);

  bool useBB = false;

  if (q != 1.0 || t != 0.0) {
    useBB = true;
  }

  if (useBB) {
    recBB(source, result, index, dim, 0.0, 0.0, q, t);
  } else {
    rec(source, result, index, dim, 0.0, 0.0);
  }
}

void XPhidPhiDownBBLinear::rec(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                               grid_iterator& index, size_t dim, double fl, double fr) {
  size_t seq = index.seq();

  double alpha_value = source[seq];

  sgpp::base::level_t l;
  sgpp::base::index_t i;

  index.get(dim, l, i);

  double hhalf = 1.0 / static_cast<double>(1 << (l + 1));
  double i_dbl = static_cast<double>(i);

  // integration
  result[seq] =
      (((fl * ((hhalf * i_dbl) - hhalf)) + (fr * (((-1.0) * (hhalf * i_dbl)) - hhalf))) -
       ((1.0 / 3.0) * (((1.0 / static_cast<double>(1 << l))) * alpha_value)));  // diagonal entry

  // dehierarchisation
  double fm = ((fl + fr) / 2.0) + alpha_value;

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

void XPhidPhiDownBBLinear::recBB(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                 grid_iterator& index, size_t dim, double fl, double fr,
                                 double q, double t) {
  size_t seq = index.seq();

  double alpha_value = source[seq];

  sgpp::base::level_t l;
  sgpp::base::index_t i;

  index.get(dim, l, i);

  double hhalf = 1.0 / static_cast<double>(1 << (l + 1));
  double i_dbl = static_cast<double>(i);

  // integration
  result[seq] = (((fl * ((q * ((hhalf * i_dbl) - hhalf)) + (0.5 * t))) +
                  (fr * ((q * (((-1.0) * (hhalf * i_dbl)) - hhalf)) - (0.5 * t)))) -
                 ((1.0 / 3.0) *
                  (((1.0 / static_cast<double>(1 << l)) * q) * alpha_value)));  // diagonal entry

  // dehierarchisation
  double fm = ((fl + fr) / 2.0) + alpha_value;

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB(source, result, index, dim, fl, fm, q, t);
    }

    index.stepRight(dim);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB(source, result, index, dim, fm, fr, q, t);
    }

    index.up(dim);
  }
}

}  // namespace finance
}  // namespace sgpp
