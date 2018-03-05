// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATASETMORTONORDER_HPP
#define DATASETMORTONORDER_HPP

#include <sgpp/datadriven/tools/Dataset.hpp>

#include <sgpp/globaldef.hpp>

#include <stdint.h>
#include <vector>

namespace sgpp {
namespace datadriven {

/// Class for re-arranging Datasets along a Morton order curve
class MortonOrder {
 public:
  /// Generates the permutation according to the given dataset.
  explicit MortonOrder(sgpp::datadriven::Dataset *dataset);

  /// Re-arrange the Dataset object along Z-Curve
  void orderDataset();
  /// Restores the original order of the Dataset object
  void restoreDataset();

  /// Check if permutation is identity
  bool isIdentity() const;

  /// Check if Dataset is orderes along Z-Curve
  bool isOrdered() const;

  /// Access to the permutation vector
  const std::vector<size_t> &getPermutation() const;

 protected:
  bool _isOrdered;
  std::vector<size_t> permutation;
  bool _isIdentity;
  sgpp::datadriven::Dataset *_dataset;
};

}  // namespace datadriven
}  // namespace sgpp

#endif  // MORTONORDER_HPP
