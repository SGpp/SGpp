// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/prewavelet/algorithm_sweep/LaplaceUpPrewavelet.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {
LaplaceUpPrewavelet::LaplaceUpPrewavelet(sgpp::base::GridStorage* storage) : storage(storage) {}

LaplaceUpPrewavelet::~LaplaceUpPrewavelet() {}

void LaplaceUpPrewavelet::operator()(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                     grid_iterator& index, size_t dim) {
  size_t seq = index.seq();
  sgpp::base::level_t l;
  sgpp::base::index_t i;
  sgpp::base::level_t l_old;
  sgpp::base::index_t i_old;
  sgpp::base::index_t last_index;
  size_t _seq;
  size_t _seql1;
  size_t _seqr1;
  size_t _seqr2;
  double _val, _vall1, _vall2, _valr1, _valr2;
  double h;
  bool hasChilds = false;
  // level_t max_level = getGridDepth(index, dim);

  index.get(dim, l, i);
  index.get(dim, l_old, i_old);

  // Level 1

  result[seq] = 1.0 / 3.0 * source[seq];

  if (!index.hintLeft(dim) && !index.hintRight(dim)) {
    return;
  }

  // Level 2
  l = 2;
  h = static_cast<double>(1 << l);
  index.set(dim, 2, 1);

  if (!hasChilds && (index.hintLeft(dim) || index.hintRight(dim))) hasChilds = true;

  _seql1 = index.seq();
  _vall1 = storage->isInvalidSequenceNumber(_seql1) ? 0.0 : source[_seql1];

  index.set(dim, 2, 3);

  if (!hasChilds && (index.hintLeft(dim) || index.hintRight(dim))) hasChilds = true;

  _seqr1 = index.seq();
  _valr1 = storage->isInvalidSequenceNumber(_seqr1) ? 0.0 : source[_seqr1];

  if (!storage->isInvalidSequenceNumber(_seql1))
    result[_seql1] = 11.0 / 75.0 * _vall1 + 1.0 / 25.0 * _valr1;

  if (!storage->isInvalidSequenceNumber(_seqr1))
    result[_seqr1] = 11.0 / 75.0 * _valr1 + 1.0 / 25.0 * _vall1;

  //
  //    result[_seql1] = 11.0 / 75.0 * source[_seql1];
  //    result[_seqr1] = 11.0 / 75.0 * source[_seqr1];

  if (!hasChilds) {
    index.set(dim, l_old, i_old);
    return;
  }

  while (true) {
    l++;
    hasChilds = false;
    h = 1.0 / (1 << l);

    last_index = (1 << (l - 1)) - 1;  // Number of Points in this level

    index.set(dim, l, 1);

    if (!hasChilds && (index.hintLeft(dim) || index.hintRight(dim))) hasChilds = true;

    _seq = index.seq();
    _val = storage->isInvalidSequenceNumber(_seq) ? 0.0 : source[_seq];

    index.set(dim, l, 3);

    if (!hasChilds && (index.hintLeft(dim) || index.hintRight(dim))) hasChilds = true;

    _seqr1 = index.seq();
    _valr1 = storage->isInvalidSequenceNumber(_seqr1) ? 0.0 : source[_seqr1];

    index.set(dim, l, 5);

    if (!hasChilds && (index.hintLeft(dim) || index.hintRight(dim))) hasChilds = true;

    _seqr2 = index.seq();
    _valr2 = storage->isInvalidSequenceNumber(_seqr2) ? 0.0 : source[_seqr2];

    if (!storage->isInvalidSequenceNumber(_seq))
      result[_seq] = 44.0 / 75.0 * h * _val      //
                     + 11.0 / 75.0 * h * _valr1  //
                     - 1.0 / 75.0 * h * _valr2;  //

    _seql1 = _seq;
    _seq = _seqr1;
    _seqr1 = _seqr2;
    _vall1 = _val;
    _val = _valr1;
    _valr1 = _valr2;
    index.set(dim, l, 7);

    if (!hasChilds && (index.hintLeft(dim) || index.hintRight(dim))) hasChilds = true;

    _seqr2 = index.seq();
    _valr2 = storage->isInvalidSequenceNumber(_seqr2) ? 0.0 : source[_seqr2];

    if (!storage->isInvalidSequenceNumber(_seq))
      result[_seq] = 11.0 / 75.0 * h * _vall1    //
                     + 18.0 / 25.0 * h * _val    //
                     + 2.0 / 15.0 * h * _valr1   //
                     - 1.0 / 75.0 * h * _valr2;  //

    // Main loop--------------------------------------

    for (i = 2; i < last_index - 1; i++) {
      _seql1 = _seq;
      _seq = _seqr1;
      _seqr1 = _seqr2;
      _vall2 = _vall1;
      _vall1 = _val;
      _val = _valr1;
      _valr1 = _valr2;
      index.set(dim, l, i * 2 + 5);

      if (!hasChilds && (index.hintLeft(dim) || index.hintRight(dim))) hasChilds = true;

      _seqr2 = index.seq();
      _valr2 = storage->isInvalidSequenceNumber(_seqr2) ? 0.0 : source[_seqr2];

      if (!storage->isInvalidSequenceNumber(_seq))
        result[_seq] = -1.0 / 75.0 * h * _vall2    //
                       + 2.0 / 15.0 * h * _vall1   //
                       + 18.0 / 25.0 * h * _val    //
                       + 2.0 / 15.0 * h * _valr1   //
                       - 1.0 / 75.0 * h * _valr2;  //
    }

    // Main loop--------------------------------------

    _seql1 = _seq;
    _seq = _seqr1;
    _seqr1 = _seqr2;
    _vall2 = _vall1;
    _vall1 = _val;
    _val = _valr1;
    _valr1 = _valr2;

    if (!storage->isInvalidSequenceNumber(_seq))
      result[_seq] = 11.0 / 75.0 * h * _valr1    //
                     + 18.0 / 25.0 * h * _val    //
                     + 2.0 / 15.0 * h * _vall1   //
                     - 1.0 / 75.0 * h * _vall2;  //

    _seql1 = _seq;
    _seq = _seqr1;
    _vall2 = _vall1;
    _vall1 = _val;
    _val = _valr1;

    if (!storage->isInvalidSequenceNumber(_seq))
      result[_seq] = 44.0 / 75.0 * h * _val      //
                     + 11.0 / 75.0 * h * _vall1  //
                     - 1.0 / 75.0 * h * _vall2;  //

    if (!hasChilds) {
      index.set(dim, l_old, i_old);
      return;
    }
  }
}
}  // namespace pde
}  // namespace sgpp
