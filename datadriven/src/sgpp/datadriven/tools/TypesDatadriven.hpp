// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TYPESDATADRIVEN_HPP
#define TYPESDATADRIVEN_HPP

#include <sgpp/globaldef.hpp>

#include <cstddef>


namespace sgpp {

namespace datadriven {

/**
 * struct to encapsulate the classsifiers quality by its
 * characteristic numbers
 */
struct ClassificatorQuality {
  /// number of true positive classified instances
  size_t truePositive_;
  /// number of true negative classified instances
  size_t trueNegative_;
  /// number of false positive classified instances
  size_t falsePositive_;
  /// number of false negative classified instances
  size_t falseNegative_;
};

/**
 * strcut to encapsulate the learner's timings
 * during training
 */
struct LearnerTiming {
  /// complete learning time
  double timeComplete_;
  /// time to apply B (including data transfer to eventually used accelerators)
  double timeMultComplete_;
  /// pure application time of B
  double timeMultCompute_;
  /// time to apply B^T (including data transfer to eventually used accelerators)
  double timeMultTransComplete_;
  /// pure application time of B^T
  double timeMultTransCompute_;
  /// time of regularization
  double timeRegularization_;
  /// number of executed Floating Point operations
  double GFlop_;
  /// number of transferred Gbytes
  double GByte_;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* TYPESDATADRIVEN_HPP */

