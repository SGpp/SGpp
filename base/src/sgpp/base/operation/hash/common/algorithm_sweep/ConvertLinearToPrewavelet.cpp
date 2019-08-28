// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/ConvertLinearToPrewavelet.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

void ConvertLinearToPrewavelet::operator()(DataVector& source,
    DataVector& result,
    grid_iterator& index, size_t dim) {

  // The Algorithm runs bottom-up, thus we need the maximal grid-depth
  // in the current dimension
  level_type max_level = index.getGridDepth(dim);

  if (max_level == 1) {
    return;
  }

  level_type level = max_level;

  level_type init_level;
  index_type init_index;

  index.get(dim, init_level, init_index);

  size_t _seq;
  double _val;

  double* temp = new double[1 << (level + 1)];  // The temp values
  double* r = new double[1 << (level - 1)];
  // The following arrays are required for the triangulation
  double* gam = new double[1 << (level - 1)];
  double* u = new double[1 << (level - 1)];

  for (int i = 0; i < 1 << (level - 1); ++i) {
    r[i] = 0.0;
    gam[i] = 0.0;
    u[i] = 0.0;
  }

  for (int i = 0; i < 1 << (level + 1); ++i) {
    temp[i] = 0.0;
  }

  for (level = max_level; level >= 2; --level) {
    // First, we set the right-hand side of the triangular eqaution system

    // special treatment for the first point in the row
    index.set(dim, level, 1);
    _seq = index.seq();
    _val = storage.isInvalidSequenceNumber(_seq) ? 0.0 : result[_seq];
    r[0] = _val - temp[1] + 0.5 * temp[3];

    // special treatment for the last point in the row
    index.set(dim, level, (1 << level) - 1);
    _seq = index.seq();
    _val = storage.isInvalidSequenceNumber(_seq) ? 0.0 : result[_seq];
    r[(1 << (level - 1)) - 1] = _val - temp[(1 << (level + 1)) - 3]
                                + 0.5 * temp[(1 << (level + 1)) - 5];

    // normal treatment for the middle
    if (level > 2) {
      for (index_type ind = 3; ind
           < static_cast<unsigned int>((1 << level) - 1); ind = ind + 2) {
        index.set(dim, level, ind);
        _seq = index.seq();
        _val = storage.isInvalidSequenceNumber(_seq) ? 0.0 : result[_seq];
        r[ind / 2] = _val - temp[2 * ind - 1] + 0.5 * temp[2 * ind
                     - 3] + 0.5 * temp[2 * ind + 1];
      }
    }

    // Run the actual triangulation

    double bet = 0.0;

    // This is the forward-reduction
    for (int i = 0; i < 1 << (level - 1); i++) {
      index.set(dim, level, i * 2 + 1);

      if (storage.isInvalidSequenceNumber(index.seq())) {
        u[i] = 0;
        gam[i] = 0;
        continue;
      } else {
        if (i == 0) {
          bet = 1.2;
          u[0] = (r[0] / bet);
          continue;
        }

        /* We need to distinguish between the start of the triangulation
         * and the normal start. Please note, this is done in that strange
         * because in the adaptive case it is perfectly possible to have
         * grid points just in the middle of the level, thus the beginning
         * of the triangulation is not determined with index 1.
         */

        index.set(dim, level, i * 2 - 1);

        if (!storage.isInvalidSequenceNumber(index.seq()))
          gam[i] = 0.4 / bet;
        else
          gam[i] = 0.0;

        if (i == (1 << (level - 1)) - 1) {
          bet = 1.2;
        } else {
          bet = 1.6;
        }

        index.set(dim, level, i * 2 + 1);
        bet = bet - 0.4 * gam[i];
        u[i] = (r[i] - 0.4 * u[i - 1]) / bet;
      }
    }

    // Backward-Reduction
    for (int i = (1 << (level - 1)) - 2; i >= 0; --i) {
      u[i] = u[i] - gam[i + 1] * u[i + 1];
    }

    // Now the results are all ready in u.

    for (int i = 0; i < 1 << (level - 1); i++) {
      index.set(dim, level, (i * 2 + 1));
      _seq = index.seq();

      if (!storage.isInvalidSequenceNumber(_seq))
        result[_seq] = u[i];
    }

    // create new temp valuesand the normal start.
    // Please note, this is done in that strange
    for (index_type ind = 1; ind < static_cast<unsigned int>((1 << level) - 2); ind
         = ind + 2) {
      index.set(dim, level, ind);
      _seq = index.seq();
      _val = storage.isInvalidSequenceNumber(_seq) ? 0.0 : result[_seq];

      if (level != max_level) {
        temp[ind] = temp[ind * 2 + 1];
      }

      temp[ind] = temp[ind] - 0.6 * _val;

      index.set(dim, level, ind + 2);
      _seq = index.seq();
      _val = storage.isInvalidSequenceNumber(_seq) ? 0.0 : result[_seq];

      temp[ind] = temp[ind] - 0.6 * _val;
    }
  }

  // Treatment of the top-point in this dimension
  index.set(dim, init_level, init_index);
  _seq = index.seq();
  _val = storage.isInvalidSequenceNumber(_seq) ? 0.0 : result[_seq];

  if (!storage.isInvalidSequenceNumber(_seq))
    result[_seq] = _val - temp[1];

  delete[] temp;
  delete[] r;
  delete[] gam;
  delete[] u;
}

}  // namespace base
}  // namespace sgpp
