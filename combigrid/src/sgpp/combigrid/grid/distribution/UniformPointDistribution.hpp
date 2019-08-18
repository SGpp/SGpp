// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/grid/distribution/AbstractPointDistribution.hpp>

namespace sgpp {
namespace combigrid {

/**
 * provides uniform points, i. e. {(k + 1)/(n + 1) for k = 0, ..., n-1}
 */
class UniformPointDistribution : public AbstractPointDistribution {
 public:
  virtual ~UniformPointDistribution();

  virtual double compute(size_t numPoints, size_t j);
};

} /* namespace combigrid */
} /* namespace sgpp*/
