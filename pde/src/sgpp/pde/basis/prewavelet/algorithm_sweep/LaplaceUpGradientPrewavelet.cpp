// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/prewavelet/algorithm_sweep/LaplaceUpGradientPrewavelet.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

LaplaceUpGradientPrewavelet::LaplaceUpGradientPrewavelet(sgpp::base::GridStorage* storage)
    : storage(storage) {}

LaplaceUpGradientPrewavelet::~LaplaceUpGradientPrewavelet() {}

void LaplaceUpGradientPrewavelet::operator()(sgpp::base::DataVector& source,
                                             sgpp::base::DataVector& result, grid_iterator& index,
                                             size_t dim) {
  sgpp::base::level_t l = index.getGridDepth(dim);
  sgpp::base::index_t i;

  sgpp::base::level_t l_old;
  sgpp::base::index_t i_old;

  index.get(dim, l_old, i_old);

  size_t _seq;
  size_t _seql;
  size_t _seqr;

  double _vall, _valr;

  double h;

  if (l == 1) return;

  index.set(dim, l, 1);
  _seqr = index.seq();
  _valr = storage->isInvalidSequenceNumber(_seqr) ? 0.0 : source[_seqr];

  double* temp_current = new double[(1 << (l - 1)) - 1];

  for (i = 0; i < static_cast<unsigned int>(1 << (l - 1)) - 1; i++) {
    _seql = _seqr;
    _vall = _valr;
    index.set(dim, l, 2 * i + 3);
    _seqr = index.seq();
    _valr = storage->isInvalidSequenceNumber(_seqr) ? 0.0 : source[_seqr];
    temp_current[i] = -0.6 * (_vall + _valr);
  }

  l--;

  for (; l > 2; l--) {
    i = 0;
    h = static_cast<double>(1 << l);
    index.set(dim, l, i + 1);
    _seq = index.seq();

    if (!storage->isInvalidSequenceNumber(_seq))
      result[_seq] =
          0.9 * h * (2 * temp_current[i] - temp_current[i + 1]) -
          0.6 * h * (-temp_current[i + 2] + 2 * temp_current[i + 1] - temp_current[i]) +
          0.1 * h * (-temp_current[i + 1] + 2 * temp_current[i + 2] - temp_current[i + 3]);

    i = 2;
    index.set(dim, l, i + 1);
    _seq = index.seq();

    if (!storage->isInvalidSequenceNumber(_seq))
      result[_seq] =
          h * (2 * temp_current[i] - temp_current[i + 1] - temp_current[i - 1]) +
          -0.6 * h * (-temp_current[i + 2] + 2 * temp_current[i + 1] - 2 * temp_current[i]) -
          0.6 * h * (-temp_current[i - 2] + 2 * temp_current[i - 1]) +
          0.1 * h * (-temp_current[i + 1] + 2 * temp_current[i + 2] - temp_current[i + 3]) +
          0.1 * h * (-temp_current[i - 1] + 2 * temp_current[i - 2]);

    for (i = 4; i < static_cast<unsigned int>(1 << l) - 4; i = i + 2) {
      index.set(dim, l, i + 1);
      _seq = index.seq();

      if (!storage->isInvalidSequenceNumber(_seq))
        result[_seq] =
            h * (2 * temp_current[i] - temp_current[i + 1] - temp_current[i - 1]) -
            0.6 * h * (-temp_current[i + 2] + 2 * temp_current[i + 1] - 2 * temp_current[i]) -
            0.6 * h * (-temp_current[i - 2] + 2 * temp_current[i - 1]) +
            0.1 * h * (-temp_current[i + 1] + 2 * temp_current[i + 2] - temp_current[i + 3]) +
            0.1 * h * (-temp_current[i - 1] + 2 * temp_current[i - 2] - temp_current[i - 3]);
    }

    i = (1 << l) - 4;
    index.set(dim, l, i + 1);
    _seq = index.seq();

    if (!storage->isInvalidSequenceNumber(_seq))
      result[_seq] =
          h * (2 * temp_current[i] - temp_current[i + 1] - temp_current[i - 1]) +
          -0.6 * h * (-temp_current[i + 2] + 2 * temp_current[i + 1] - 2 * temp_current[i]) -
          0.6 * h * (-temp_current[i - 2] + 2 * temp_current[i - 1]) +
          0.1 * h * (-temp_current[i + 1] + 2 * temp_current[i + 2]) +
          0.1 * h * (-temp_current[i - 1] + 2 * temp_current[i - 2] - temp_current[i - 3]);

    i = (1 << l) - 2;
    index.set(dim, l, i + 1);
    _seq = index.seq();

    if (!storage->isInvalidSequenceNumber(_seq))
      result[_seq] =
          0.9 * h * (2 * temp_current[i] - temp_current[i - 1]) -
          0.6 * h * (-temp_current[i - 2] + 2 * temp_current[i - 1] - temp_current[i]) +
          0.1 * h * (-temp_current[i - 1] + 2 * temp_current[i - 2] - temp_current[i - 3]);

    index.set(dim, l, 1);
    _seqr = index.seq();
    _valr = storage->isInvalidSequenceNumber(_seqr) ? 0.0 : source[_seqr];

    for (i = 0; i < static_cast<unsigned int>(1 << (l - 1)) - 1; i++) {
      _seql = _seqr;
      _vall = _valr;
      index.set(dim, l, 2 * i + 3);
      _seqr = index.seq();
      _valr = storage->isInvalidSequenceNumber(_seqr) ? 0.0 : source[_seqr];
      temp_current[i] = -0.6 * (_vall + _valr) + temp_current[2 * i + 1];
    }
  }

  if (l == 2) {
    h = static_cast<double>(1 << l);

    index.set(dim, 2, 1);
    _seql = index.seq();

    if (!storage->isInvalidSequenceNumber(_seql)) {
      result[_seql] = 0.9 * h * (2 * temp_current[0] - temp_current[1]) -
                      0.6 * h * (-temp_current[2] + 2 * temp_current[1] - temp_current[0]) +
                      0.1 * h * (-temp_current[1] + 2 * temp_current[2]);
      _vall = source[_seql];
    } else {
      _vall = 0.0;
    }

    index.set(dim, 2, 3);
    _seqr = index.seq();

    if (!storage->isInvalidSequenceNumber(_seqr)) {
      result[_seqr] = 0.9 * h * (2 * temp_current[2] - temp_current[1]) -
                      0.6 * h * (-temp_current[0] + 2 * temp_current[1] - temp_current[2]) +
                      0.1 * h * (-temp_current[1] + 2 * temp_current[0]);
      _valr = source[_seqr];
    } else {
      _valr = 0.0;
    }

    temp_current[0] = -0.6 * (_vall + _valr) + temp_current[1];

    l = 1;
  }

  if (l == 1) {
    h = static_cast<double>(1 << l);
    index.set(dim, 1, 1);
    _seq = index.seq();
    result[_seq] = 2 * h * temp_current[0];
  }

  // I dont think thats necessary, but just to be sure!
  index.set(dim, l_old, i_old);
  delete[] temp_current;
}
}  // namespace pde
}  // namespace sgpp
