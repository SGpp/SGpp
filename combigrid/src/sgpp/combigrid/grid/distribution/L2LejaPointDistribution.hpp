// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/distribution/AbstractPointDistribution.hpp>

#include <cmath>
#include <functional>
#include <vector>

namespace sgpp {
namespace combigrid {

class L2LejaPointDistribution : public AbstractPointDistribution {
  std::vector<double> points;
  std::vector<double> sortedPoints;  // also include boundaries 0.0 and 1.0
  SingleFunction weightFunction;
  size_t numAdditionalPoints;

  void addPoint(double point);

  void computeNextPoint();
  void computeNextPointLogSumExp();

 public:
  L2LejaPointDistribution();
  explicit L2LejaPointDistribution(SingleFunction weightFunction, size_t numAdditionalPoints = 10);
  virtual ~L2LejaPointDistribution();

  virtual double compute(size_t numPoints, size_t j);
};

} /* namespace combigrid */
} /* namespace sgpp */
