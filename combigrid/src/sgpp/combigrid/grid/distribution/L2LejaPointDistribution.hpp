// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/distribution/AbstractPointDistribution.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>

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

  /**
   * This method computes the next L2-Leja point given the current sequence. This method uses
   * the log-sum-exp trick to compute the relevant integrals. The reason is that the integrals
   * become very small such that they do not fit anoymore into a double. This tricks avoids this
   * underflow issue.
   */
  void computeNextPoint();

 public:
  L2LejaPointDistribution();
  explicit L2LejaPointDistribution(SingleFunction weightFunction, size_t numAdditionalPoints = 10);
  virtual ~L2LejaPointDistribution();

  virtual double compute(size_t numPoints, size_t j);

 private:
  double evaluate_denominator(size_t& argmaxIndex, GaussLegendreQuadrature& quadRule,
                              double tol = 1e-14);
  double evaluate_numerator(size_t argmaxIndex, GaussLegendreQuadrature& quadRule,
                            double tol = 1e-14);
};

} /* namespace combigrid */
} /* namespace sgpp */
