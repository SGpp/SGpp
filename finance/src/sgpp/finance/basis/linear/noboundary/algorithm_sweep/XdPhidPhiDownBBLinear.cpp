// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XdPhidPhiDownBBLinear.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace finance {



XdPhidPhiDownBBLinear::XdPhidPhiDownBBLinear(SGPP::base::GridStorage* storage) :
  storage(storage), boundingBox(storage->getBoundingBox()) {
}

XdPhidPhiDownBBLinear::~XdPhidPhiDownBBLinear() {
}

void XdPhidPhiDownBBLinear::operator()(SGPP::base::DataVector& source,
                                       SGPP::base::DataVector& result, grid_iterator& index, size_t dim) {
  float_t q = boundingBox->getIntervalWidth(dim);
  float_t t = boundingBox->getIntervalOffset(dim);

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

void XdPhidPhiDownBBLinear::rec(SGPP::base::DataVector& source,
                                SGPP::base::DataVector& result, grid_iterator& index, size_t dim, float_t fl,
                                float_t fr) {
  size_t seq = index.seq();

  float_t alpha_value = source[seq];

  SGPP::base::GridStorage::index_type::level_type l;
  SGPP::base::GridStorage::index_type::index_type i;

  index.get(dim, l, i);

  float_t i_dbl = static_cast<float_t>(i);
  //int l_int = static_cast<int>(l);

  //float_t h = (1.0/(static_cast<float_t>(1<<(l_int))));

  float_t diagonal = 2 * i_dbl;

  // integration
  result[seq] = (  (0.5 * (fl - fr)) + (diagonal * alpha_value) );

  // dehierarchisation
  float_t fm = ((fl + fr) / 2.0) + alpha_value;

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

void XdPhidPhiDownBBLinear::recBB(SGPP::base::DataVector& source,
                                  SGPP::base::DataVector& result, grid_iterator& index, size_t dim, float_t fl,
                                  float_t fr, float_t q, float_t t) {
  size_t seq = index.seq();

  float_t alpha_value = source[seq];

  SGPP::base::GridStorage::index_type::level_type l;
  SGPP::base::GridStorage::index_type::index_type i;

  index.get(dim, l, i);

  //float_t i_dbl = static_cast<float_t>(i);
  int l_int = static_cast<int>(l);

  float_t h = (1.0 / (static_cast<float_t>(1 << (l_int))));

  float_t diagonal = (2 * (static_cast<float_t>(i) * h * q + t)) / (q * h);

  // integration
  result[seq] = (  (0.5 * q * (fl - fr)) + (diagonal * alpha_value) );

  // dehierarchisation
  float_t fm = ((fl + fr) / 2.0) + alpha_value;

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->end(index.seq())) {
      recBB(source, result, index, dim, fl, fm, q, t);
    }

    index.stepRight(dim);

    if (!storage->end(index.seq())) {
      recBB(source, result, index, dim, fm, fr, q, t);
    }

    index.up(dim);
  }
}

// namespace detail

} // namespace SGPP
}
