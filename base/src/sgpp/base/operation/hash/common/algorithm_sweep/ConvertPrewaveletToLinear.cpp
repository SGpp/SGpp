// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/ConvertPrewaveletToLinear.hpp>


#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

ConvertPrewaveletToLinear::ConvertPrewaveletToLinear(GridStorage& storage) :
  storage(storage) {
}

ConvertPrewaveletToLinear::~ConvertPrewaveletToLinear() {
}

void ConvertPrewaveletToLinear::operator()(DataVector& source,
    DataVector& result,
    grid_iterator& index, size_t dim) {
  level_type max_level = index.getGridDepth(dim);

  if (max_level == 1) {
    return;
  }

  level_type level = max_level;

  level_type init_level;
  index_type init_index;
  size_t _seq;
  size_t _seq_temp;
  double _val = 0.0;
  double* temp_current = nullptr;
  double* temp_old = nullptr;


  index.get(dim, init_level, init_index);

  for (; level > 1; --level) {
    if (level == max_level) {
      temp_current = new double[1 << level];
      temp_old = new double[1 << (level + 1)];

      for (int t = 0; t < (1 << (level + 1)); t++) {
        temp_old[t] = 0;
      }

      for (int t = 2; t < (1 << level); t = t + 2) {
        index.set(dim, level, t - 1);
        _seq = index.seq();
        _val = storage.isInvalidSequenceNumber(_seq) ? 0.0 : source[_seq];
        temp_current[t] = -0.6 * _val;

        index.set(dim, level, t + 1);
        _seq = index.seq();
        _val = storage.isInvalidSequenceNumber(_seq) ? 0.0 : source[_seq];
        temp_current[t] = temp_current[t] - 0.6 * _val;
      }
    } else {
      delete[] temp_old;
      temp_old = temp_current;
      temp_current = new double[(1 << level)];

      for (int t = 2; t < (1 << level); t = t + 2) {
        index.set(dim, level, t - 1);
        _seq = index.seq();
        _val = storage.isInvalidSequenceNumber(_seq) ? 0.0 : source[_seq];
        temp_current[t] = -0.6 * _val + temp_old[t * 2];

        index.set(dim, level, t + 1);
        _seq = index.seq();
        _val = storage.isInvalidSequenceNumber(_seq) ? 0.0 : source[_seq];
        temp_current[t] = temp_current[t] - 0.6 * _val;
      }
    }

    // Special treatment for first index
    index.set(dim, level, 1);
    _seq = index.seq();
    _val = storage.isInvalidSequenceNumber(_seq) ? 0.0 : source[_seq];
    double current_value = _val;
    double left_value = 0;

    if (!storage.isInvalidSequenceNumber(_seq))
      result[_seq] = 0.9 * current_value;

    index.set(dim, level, 3);
    _seq_temp = index.seq();
    _val = storage.isInvalidSequenceNumber(_seq_temp) ? 0.0 : source[_seq_temp];

    if (!storage.isInvalidSequenceNumber(_seq)) {
      result[_seq] += 0.1 * _val;
      result[_seq] += temp_old[2] - 0.5 * temp_current[2];
    }

    for (int i = 3; i < (1 << level) - 2; i = i + 2) {
      index.set(dim, level, i);
      _seq = index.seq();

      index.set(dim, level, i + 2);
      _seq_temp = index.seq();

      left_value = current_value;
      _val = storage.isInvalidSequenceNumber(_seq) ? 0.0 : source[_seq];
      current_value = _val;

      _val = storage.isInvalidSequenceNumber(_seq_temp) ? 0.0 : source[_seq_temp];

      if (!storage.isInvalidSequenceNumber(_seq)) {
        result[_seq] = current_value + 0.1 * left_value + 0.1
                       * _val;

        result[_seq] += temp_old[i * 2] - 0.5 * temp_current[i - 1]
                        - 0.5 * temp_current[i + 1];
      }
    }

    // Special treatment for last index
    index_type last = (1 << static_cast<index_type>(level)) - 1;
    index.set(dim, level, last);
    _seq = index.seq();
    _val = storage.isInvalidSequenceNumber(_seq) ? 0.0 : source[_seq];

    if (!storage.isInvalidSequenceNumber(_seq)) {
      result[_seq] = 0.9 * _val + 0.1 * current_value;

      result[_seq] += temp_old[last * 2] - 0.5 * temp_current[last
                      - 1];
    }
  }

  index.set(dim, init_level, init_index);
  _seq = index.seq();

  result[_seq] = source[_seq] + temp_current[2];

  delete[] temp_old;
  delete[] temp_current;
}

}  // namespace base
}  // namespace sgpp
