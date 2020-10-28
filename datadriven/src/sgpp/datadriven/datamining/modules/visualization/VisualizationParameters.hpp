// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include<cstddef>

namespace sgpp {
namespace datadriven {

struct VisualizationParameters {
  /**
   * The perplexity to use in case tsne is the selected algorithm
   */
  double perplexity_ = 30;

  /**
   * The theta parameter to use in case tsne is the selected algorithm
   */
  double theta_ = 0.5;

  /*
   * The random seed to initialize the selected algorithm
   */
  std::size_t seed_ = 100;

  /*
   * The maximum number of iteration to run on the gradient descent of a selected
   * algorithm
   */
  std::size_t maxNumberIterations_ = 1000;

  /*
   * The dimentionality to which we want to reduce the data for visualization
   * purposes
   */
  std::size_t targetDimension_ = 2;

  /*
   * Number of cores to run the algorithm in parallel
   */
  std::size_t numberCores_ = 1;
};
}  // namespace datadriven
}  // namespace sgpp
