// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SORTEDDATASET_HPP
#define SORTEDDATASET_HPP

#include <sgpp/datadriven/tools/Dataset.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

/// Dataset that can be ordered. Accessing the included DataMatrix might invalidate the order.
class SortedDataset : public Dataset {
 public:
  /// Available permutations
  enum OrderType { None, External, Random, Morton, Invalid };

  /**
   * Constructs an empty dataset (zero size).
   */
  SortedDataset();

  /**
   * Constructs an empty dataset with given size.
   *
   * @param numberInstances number of instances in the dataset
   * @param dimension number of dimensions in the dataset
   */
  SortedDataset(size_t numberInstances, size_t dimension);

  /**
   * Constructs a copy of a dataset.
   */
  explicit SortedDataset(const Dataset &src);

  /** Sets the OrderType to OrderType::Invalid.
   *  @return training data of the dataset
   */
  sgpp::base::DataMatrix &getData();

  /** Sets the order for the dataset and rearranges the data.
   *  The order is unchanged in case of OrderType::External.
   *  OrderType::None creates the identity permutation.
   */
  void setOrder(OrderType order);

  /** Reorders the data with a given permutation.
   *  Sets the order type to OrderType::External
   *  If the permutation is invalid due to a size or index mismatch, ot is set to
   * OrderType::Invalid.
   *  The i-th data value is permuted to the permutation[i]-th position.
   */
  void setOrder(const std::vector<size_t> &permutation);

  /// Returns current order type
  OrderType getOrderType() const;

  /// Restores the original order in case of a valid order type. Sets ot to OrderType::None
  void restoreOrder();

 protected:
  OrderType ot;
  std::vector<size_t> perm;

  /// uses the permutation on the DataMatrix and the DataVector
  void usePermutation();
};

}  // namespace datadriven
}  // namespace sgpp

#endif
