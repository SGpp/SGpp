// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/tools/SortedDataset.hpp>

#include <algorithm>  // std::random_shuffle
#include <cstdlib>    // std::rand, std::srand
#include <ctime>      // std::time
#include <vector>     // std::vector

namespace sgpp {
namespace datadriven {

///@cond DOXY_IGNORE // NOLINT()
namespace SortedDatasetDetail {

union ext_double_t {
  double val;  // Value
  struct {
    uint64_t dig : 52;  // Digits
    uint32_t exp : 11;  // Exponent
    uint8_t sig : 1;    // Sign
  } bit;
};

/// Returns most significant bit from a integer value
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

/// Returns MSB from 2 double values
int XOR_MSB(ext_double_t a, ext_double_t b) {
  int ret;
  if (a.bit.exp == b.bit.exp) {
    ret = a.bit.exp + MSB(a.bit.dig ^ b.bit.dig);
  } else if (a.bit.exp > b.bit.exp) {
    ret = a.bit.exp + 52;
  } else {
    ret = b.bit.exp + 52;
  }
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

}  // namespace SortedDatasetDetail
///@endcond // NOLINT()

using SortedDatasetDetail::zorder;

/**
 * Constructs an empty dataset (zero size).
 */
SortedDataset::SortedDataset() : Dataset() {
  ot = OrderType::None;
  perm.clear();
  std::srand(static_cast<unsigned>(std::time(0)));
}

/**
 * Constructs an empty dataset with given size.
 *
 * @param numberInstances number of instances in the dataset
 * @param dimension number of dimensions in the dataset
 */
SortedDataset::SortedDataset(size_t numberInstances, size_t dimension)
    : Dataset(numberInstances, dimension) {
  ot = OrderType::None;
  perm.resize(numberInstances);
  for (size_t i = 0; i < perm.size(); ++i) perm[i] = i;
  std::srand(static_cast<unsigned>(std::time(0)));
}

/**
 * Constructs a copy of a dataset.
 */
SortedDataset::SortedDataset(const Dataset &src) : Dataset(src) {
  ot = OrderType::None;
  perm.resize(src.getNumberInstances());
  for (size_t i = 0; i < perm.size(); ++i) perm[i] = i;
  std::srand(static_cast<unsigned>(std::time(0)));
}

/** Sets the OrderType to OrderType::Invalid.
 *  @return training data of the dataset
 */
sgpp::base::DataMatrix &SortedDataset::getData() {
  ot = OrderType::Invalid;
  return Dataset::getData();
}

/** Sets the order for the dataset and rearranges the data.
 *  The order is unchanged in case of OrderType::External.
 */
void SortedDataset::setOrder(OrderType order) {
  switch (order) {
    case OrderType::None:
      ot = OrderType::None;
      perm.resize(numberInstances);
      for (size_t i = 0; i < perm.size(); ++i) perm[i] = i;
      break;
    case OrderType::Random:
      ot = OrderType::Random;
      perm.resize(numberInstances);
      for (size_t i = 0; i < perm.size(); ++i) perm[i] = i;
      std::random_shuffle(perm.begin(), perm.end());
      usePermutation();
      break;
    case OrderType::Morton:
      ot = OrderType::Morton;
      zorder(data, perm);
      usePermutation();
      break;
    case OrderType::External:
    case OrderType::Invalid:
      break;
  }
}

/** Reorders the data with a given permutation.
 *  Sets the order type to OrderType::External
 *  If the permutation is invalid due to a size or index mismatch, ot is set to OrderType::Invalid.
 *  The i-th data value is permuted to the permutation[i]-th position.
 */
void SortedDataset::setOrder(const std::vector<size_t> &permutation) {
  std::vector<size_t> check(numberInstances, 0);
  bool isValid = true;
  if (permutation.size() != numberInstances) isValid = false;
  if (isValid) {
    for (size_t i = 0; i < numberInstances; ++i) {
      if (permutation[i] > numberInstances) {
        isValid = false;
        break;
      }
      if ((++check[permutation[i]]) > 1) {
        isValid = false;
        break;
      }
    }
  }
  if (isValid) {
    ot = OrderType::External;
    perm = permutation;
    usePermutation();
  } else {
    ot = OrderType::Invalid;
  }
}

/// Returns current order type
SortedDataset::OrderType SortedDataset::getOrderType() const { return ot; }

/// Restores the original order in case of a valid order type. Sets ot to OrderType::None
void SortedDataset::restoreOrder() {
  if (ot == OrderType::Invalid) return;
  sgpp::base::DataMatrix matrix(data);
  sgpp::base::DataVector vector(targets);
  for (size_t i = 0; i < perm.size(); ++i) {
    targets[perm[i]] = vector[i];
    for (size_t d = 0; d < matrix.getNcols(); ++d) {
      data(perm[i], d) = matrix(i, d);
    }
  }
}

/// uses the permutation on the DataMatrix and the DataVector
void SortedDataset::usePermutation() {
  base::DataMatrix matrix(data);
  base::DataVector vector(targets);
  for (size_t i = 0; i < perm.size(); ++i) {
    targets[i] = vector[perm[i]];
    for (size_t d = 0; d < matrix.getNcols(); ++d) {
      data(i, d) = matrix(perm[i], d);
    }
  }
}

}  // namespace datadriven
}  // namespace sgpp
