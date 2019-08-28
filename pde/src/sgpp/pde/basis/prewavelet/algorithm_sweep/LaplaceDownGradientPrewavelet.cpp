// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/prewavelet/algorithm_sweep/LaplaceDownGradientPrewavelet.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

LaplaceDownGradientPrewavelet::LaplaceDownGradientPrewavelet(sgpp::base::GridStorage* storage)
    : storage(storage) {}

LaplaceDownGradientPrewavelet::~LaplaceDownGradientPrewavelet() {}

void LaplaceDownGradientPrewavelet::operator()(sgpp::base::DataVector& source,
                                               sgpp::base::DataVector& result, grid_iterator& index,
                                               size_t dim) {
  size_t _seql2, _seql1, _seq, _seqr1, _seqr2 = 0;
  double _vall2, _vall1, _val, _valr1, _valr2 = 0;

  sgpp::base::level_t l, l_old, max_level;
  sgpp::base::index_t i, i_old;

  double h;

  index.get(dim, l_old, i_old);         // save iterator and restore afterwards
  max_level = index.getGridDepth(dim);  // How many levels to decent?

  // We need the temp_values of the level above
  double* temp_current = new double[1];
  double* temp_old;

  // Level 1-----------------------------------------------------------
  _seq = index.seq();
  temp_current[0] = 4 * source[_seq];
  result[_seq] = 4 * source[_seq];

  // Is there a second level?
  if (max_level == 1) {
    // No, delete temp and return
    delete[] temp_current;
    index.set(dim, l_old, i_old);
    return;
  }

  // Level 2-----------------------------------------------------------
  l = 2;
  h = static_cast<double>(1 << l);
  temp_old = temp_current;
  temp_current = new double[3];
  index.leftChild(dim);
  _seql1 = index.seq();
  index.stepRight(dim);
  _seqr1 = index.seq();

  // Calculate new temp variables
  temp_current[0] = 2.4 * h * source[_seql1] + 0.8 * h * source[_seqr1];
  temp_current[1] = temp_old[0] - 2.2 * h * (source[_seql1] + source[_seqr1]);
  temp_current[2] = 2.4 * h * source[_seqr1] + 0.8 * h * source[_seql1];

  // Calculate results
  _vall1 = source[_seql1];  // We save the variables to enable work in sito.
  _valr1 = source[_seqr1];
  result[_seql1] = -0.6 * temp_old[0] + 3.56 * h * _vall1 + 2.28 * h * _valr1;
  result[_seqr1] = -0.6 * temp_old[0] + 3.56 * h * _valr1 + 2.28 * h * _vall1;

  // All other levels--------------------------------------------------
  for (l = 3; l <= max_level; l++) {
    h = static_cast<double>(1 << l);
    index.set(dim, l, 1);

    // Calculate new temp-Variables ################################################
    delete[] temp_old;
    temp_old = temp_current;
    temp_current = new double[(1 << l) - 1];

    if (l != max_level) {  // If we are in the last iteration, we don't need any more
                           // temp-variables, thus skip the calculation
      // Ramp up
      _seql2 = index.seq();
      _vall2 = storage->isInvalidSequenceNumber(_seql2) ? 0.0 : source[_seql2];

      index.stepRight(dim);
      _seql1 = index.seq();
      _vall1 = storage->isInvalidSequenceNumber(_seql1) ? 0.0 : source[_seql1];

      index.stepRight(dim);
      _seqr1 = index.seq();
      _valr1 = storage->isInvalidSequenceNumber(_seqr1) ? 0.0 : source[_seqr1];

      index.stepRight(dim);
      _seqr2 = index.seq();
      _valr2 = storage->isInvalidSequenceNumber(_seqr2) ? 0.0 : source[_seqr2];

      temp_current[0] = 2.4 * h * _vall2     // sgpp::base::Grid-point
                        + 0.8 * h * _vall1;  // right neighbor

      temp_current[1] = temp_old[0]          //
                        - 2.2 * h * _vall2   // left neighbor
                        - 2.3 * h * _vall1   // right neighbor
                        - 0.1 * h * _valr1;  // right-right neighbor

      for (i = 2; static_cast<int>(i) < (1 << l) - 3; i++) {
        if (i % 2 == 0) {  // On sgpp::base::Grid-Points
          if (static_cast<int>(i) == (1 << l) - 4) {
            temp_current[i] = 3.2 * h * _valr1                // current grid point
                              + 0.8 * h * (_vall1 + _valr2);  // neighbors left and right
          } else {
            temp_current[i] = 3.2 * h * _vall1                // current grid point
                              + 0.8 * h * (_vall2 + _valr1);  // neighbors left and right
          }
        } else {                                            // between grid points
          temp_current[i] = temp_old[(i - 1) / 2]           // temp-value from above
                            - 2.3 * h * (_vall1 + _valr1)   // left and right neighbor
                            - 0.1 * h * (_vall2 + _valr2);  // left-left and right-right neighbor

          // forward indoces except in the last iteration
          if (static_cast<int>(i) != (1 << l) - 5) {
            _vall2 = _vall1;
            _vall1 = _valr1;
            _valr1 = _valr2;
            index.stepRight(dim);
            _seq = index.seq();
            _valr2 = storage->isInvalidSequenceNumber(_seq) ? 0.0 : source[_seq];
          }
        }
      }

      // Two temp-values left
      temp_current[(1 << l) - 3] = temp_old[((1 << l) - 4) / 2]           //
                                   - 2.2 * h * _valr2 - 2.3 * h * _valr1  // left and right neighbor
                                   - 0.1 * h * _vall1;                    // left-left neighbor

      temp_current[(1 << l) - 2] = 2.4 * h * _valr2 + 0.8 * h * _valr1;  // On
                                                                         // sgpp::base::Grid-point
    }

    // End of temp-value calculation ################################################

    // Calculate result ################################################
    // Ramp-up
    index.set(dim, l, 1);
    _seql2 = index.seq();
    _vall2 = storage->isInvalidSequenceNumber(_seql2) ? 0.0 : source[_seql2];

    index.stepRight(dim);
    _seql1 = index.seq();
    _vall1 = storage->isInvalidSequenceNumber(_seql1) ? 0.0 : source[_seql1];

    index.stepRight(dim);
    _seq = index.seq();
    _val = storage->isInvalidSequenceNumber(_seq) ? 0.0 : source[_seq];

    index.stepRight(dim);
    _seqr1 = index.seq();
    _valr1 = storage->isInvalidSequenceNumber(_seqr1) ? 0.0 : source[_seqr1];

    if (l != 3) {
      index.stepRight(dim);
      _seqr2 = index.seq();
      _valr2 = storage->isInvalidSequenceNumber(_seqr2) ? 0.0 : source[_seqr2];
    }

    if (!storage->isInvalidSequenceNumber(_seql2))
      result[_seql2] = -0.6 * temp_old[0] + 3.56 * h * _vall2  // current point
                       + 2.42 * h * _vall1                     // right neighbor
                       + 0.14 * h * _val;                      // right-right neighbor

    if (!storage->isInvalidSequenceNumber(_seql1))
      result[_seql1] = -0.6 * (temp_old[0] + temp_old[1])  //
                       + 6.12 * h * _vall1                 // current point
                       + 2.42 * h * _vall2                 // left neighbor
                       + 2.56 * h * _val                   // right neighbor
                       + 0.14 * h * _valr1;                // right-right neighbor

    if (l == 3) {
      /* level = 3, that means we have only 4 gridpoints, thus we have to
       * shift the indices so that ramp-down work correctly. Otherwise
       * proceed with the inner-points
       */
      _valr2 = _valr1;
      _seqr2 = _seqr1;

      _valr1 = _val;
      _seqr1 = _seq;

      _val = _vall1;
      _seq = _seql1;

      _vall1 = _vall2;
      _seql1 = _seql2;
    } else {
      for (i = 5; static_cast<int>(i) < (1 << l) - 4; i += 2) {
        if (!storage->isInvalidSequenceNumber(_seq))
          result[_seq] = -0.6 * (temp_old[(i - 3) / 2] + temp_old[(i - 1) / 2])  //
                         + 6.12 * h * _val                                       // current point
                         + 2.56 * h * (_vall1 + _valr1)   // left and right neighbor
                         + 0.14 * h * (_vall2 + _valr2);  // left-left and right-right neighbor

        if (static_cast<int>(i) != (1 << l) - 5) {
          // forward indices, except during the last iteration
          _vall2 = _vall1;
          _vall1 = _val;
          _val = _valr1;
          _seq = _seqr1;
          _valr1 = _valr2;
          _seqr1 = _seqr2;
          index.stepRight(dim);
          _seqr2 = index.seq();
          _valr2 = storage->isInvalidSequenceNumber(_seqr2) ? 0.0 : source[_seqr2];
        }
      }
    }

    // Again, two points left
    if (!storage->isInvalidSequenceNumber(_seqr1))
      result[_seqr1] = -0.6 * (temp_old[(1 << (l - 1)) - 3] + temp_old[(1 << (l - 1)) - 2])  //
                       + 6.12 * h * _valr1   // current point
                       + 2.42 * h * _valr2   // right neighbor
                       + 2.56 * h * _val     // left neighbor
                       + 0.14 * h * _vall1;  // left-left neighbor

    if (!storage->isInvalidSequenceNumber(_seqr2))
      result[_seqr2] = -0.6 * temp_old[(1 << (l - 1)) - 2]  //
                       + 3.56 * h * _valr2                  // current point
                       + 2.42 * h * _valr1                  // left neighbor
                       + 0.14 * h * _val;                   // left-left neighbor
  }

  index.set(dim, l_old, i_old);
  delete[] temp_current;
  delete[] temp_old;
}
}  // namespace pde
}  // namespace sgpp
