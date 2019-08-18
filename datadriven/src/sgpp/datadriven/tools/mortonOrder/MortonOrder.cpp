// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/tools/mortonOrder/MortonOrder.hpp>

#include <stdint.h>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

// /@cond DOXY_IGNORE // NOLINT()
namespace MortonOrderDetail {

union ext_double_t {
  double val;  // Value
  struct {
    uint64_t dig : 52;  // Digits
    uint32_t exp : 11;  // Exponent
    uint8_t sig : 1;    // Sign
  } bit;
};

// / Returns most significant bit from a integer value
int MSB(size_t value) {
  int ret = 0;
  while (value > 1) {
    value >>= 1;
    ++ret;
  }
  return ret;
}

struct data_perm_t {
  size_t idx;
  std::vector<ext_double_t> pos;
};

// / Returns MSB from 2 double values
int XOR_MSB(ext_double_t a, ext_double_t b) {
  int ret;
  if (a.bit.exp == b.bit.exp)
    ret = a.bit.exp + MSB(a.bit.dig ^ b.bit.dig);
  else if (a.bit.exp > b.bit.exp)
    ret = a.bit.exp + 52;
  else
    ret = b.bit.exp + 52;
  return ret;
}

bool operator<(const data_perm_t &a, const data_perm_t &b) {
  size_t d = 0;
  size_t y, tmp;
  tmp = 0;
  // search the most differing dimension
  for (size_t i = 0; i < a.pos.size(); i++) {
    if ((a.pos[i].val < 0) != (b.pos[i].val < 0)) return a.pos[i].val < b.pos[i].val;

    y = XOR_MSB(a.pos[i], b.pos[i]);
    if (tmp < y) {
      tmp = y;
      d = i;
    }
  }
  // compare values in this dimension
  return a.pos[d].val < b.pos[d].val;
}

void zorder(const sgpp::base::DataMatrix &data, std::vector<size_t> &perm) {
  perm.clear();
  perm.resize(data.getNrows());
  std::vector<data_perm_t> workdata(perm.size());
  // initialize permutation as identity
  for (size_t i = 0; i < workdata.size(); ++i) {
    workdata[i].idx = i;
    workdata[i].pos.resize(data.getNcols());
    for (size_t d = 0; d < data.getNcols(); ++d) {
      workdata[i].pos[d].val = data(i, d);
    }
  }
  std::stable_sort(workdata.begin(), workdata.end());
  for (size_t i = 0; i < workdata.size(); ++i) {
    perm[i] = workdata[i].idx;
  }
}
}  // namespace MortonOrderDetail
// /@endcond  // NOLINT()

using MortonOrderDetail::zorder;

// / Constructor. Generates the permuation list on the GPU
MortonOrder::MortonOrder(sgpp::datadriven::Dataset *dataset) : _dataset(dataset) {
  // Compute permuation
  zorder(_dataset->getData(), permutation);

  // Check for identity permutation
  _isIdentity = true;
  for (size_t i = 0; i < permutation.size() && _isIdentity; ++i)
    if (permutation[i] != i) _isIdentity = false;

  _isOrdered = false;
}

// / Re-arrange a Dataset object along Z-Curve inplace
void MortonOrder::orderDataset() {
  if (_isOrdered) return;
  sgpp::base::DataMatrix matrix(_dataset->getData());
  sgpp::base::DataVector vector(_dataset->getTargets());
  for (size_t i = 0; i < permutation.size(); ++i) {
    _dataset->getTargets()[i] = vector[permutation[i]];
    for (size_t d = 0; d < matrix.getNcols(); ++d) {
      _dataset->getData()(i, d) = matrix(permutation[i], d);
    }
  }
  _isOrdered = true;
}

// / Restores the original order of a Dataset object inplace
void MortonOrder::restoreDataset() {
  if (!_isOrdered) return;
  sgpp::base::DataMatrix matrix(_dataset->getData());
  sgpp::base::DataVector vector(_dataset->getTargets());
  for (size_t i = 0; i < permutation.size(); ++i) {
    _dataset->getTargets()[permutation[i]] = vector[i];
    for (size_t d = 0; d < matrix.getNcols(); ++d) {
      _dataset->getData()(permutation[i], d) = matrix(i, d);
    }
  }
  _isOrdered = false;
}

// / Check if permutation is identity
bool MortonOrder::isIdentity() const { return _isIdentity; }

// / Check if Dataset is ordered
bool MortonOrder::isOrdered() const { return _isOrdered; }

// / Access to the permutation vector
const std::vector<size_t> &MortonOrder::getPermutation() const { return permutation; }
}  // namespace datadriven
}  // namespace sgpp
