// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/DistributionBeta.hpp>
#include <sgpp/base/tools/DistributionNormal.hpp>
#include <sgpp/base/tools/DistributionTruncExponential.hpp>
#include <sgpp/base/tools/DistributionTruncGamma.hpp>
#include <sgpp/base/tools/DistributionUniform.hpp>
#include <sgpp/base/tools/DistributionsVector.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>
#include <sgpp/base/tools/sle/solver/Eigen.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>
#include <string>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {
enum class distributionType { Normal, Beta, Uniform, Exponential, Gamma };
/**
 * a PCE object providing different methods to calculate the PCE coefficients and evaluating the
 * expansion
 */
class PolynomialChaosExpansion {
  std::function<double(const base::DataVector&)> func;
  int order;
  std::vector<distributionType> types;
  sgpp::base::DistributionsVector distributions;
  sgpp::base::DistributionsVector standardvec;
  std::vector<std::pair<double, double>> ranges;
  base::DataVector alpha;
  base::DataVector beta;
  base::DataVector coefficients;

 private:
  double evalLegendre(int n, double x);
  double evalHermite(int n, double x);
  double evalLaguerre(int n, double x);
  double evalJacobi(int n, double x, size_t i);
  double evalGenLaguerre(int n, double x, size_t i);
  std::vector<std::function<double(double, size_t)>> weights;
  std::vector<std::function<double(double, size_t)>> denoms;
  std::vector<std::function<double(double, double, size_t)>> evals;
  std::vector<std::vector<int>> multiIndex(int dimension, int order);

 public:
  /*
   * Constructor
   *
   * constructs a PCE using total-order expansion for the given function, expansion order and
   * underlying distributions
   */
  PolynomialChaosExpansion(std::function<double(const base::DataVector&)> func, int order,
                           sgpp::base::DistributionsVector distributions);

  /**
   * Destructor
   */
  ~PolynomialChaosExpansion();

  double sparseGridQuadrature(const std::function<double(const base::DataVector&)>& funct, int dim,
                              int n, size_t quadOrder);
  double adaptiveQuadratureWeighted(const std::function<double(const base::DataVector&)>& funct,
                                    int dim, size_t n, size_t quadOrder);

  /*
   * calculates the coefficients using n points and the given method
   */
  base::DataVector calculateCoefficients(int n, bool use_adaptive);
  /*
   * returns the calculated coefficients
   */
  base::DataVector getCoefficients();

  /*
   * deletes the stored coefficients
   */
  void clearCoefficients();
  /*
   * evaluates the PCE at the given point
   */
  double evalExpansion(const base::DataVector& xi, int n, bool use_adaptive);
  /*
   * returns the mean of the expansion
   */
  double getMean(int n, bool use_adaptive);
  /*
   * returns the variance of the expansion
   */
  double getVariance(int n, bool use_adaptive);
  /*
   * returns the L2 approximation error of the expansion(calculated wrt the underlying
   * distributions)
   */
  double getL2Error(int n, bool use_adaptive);
};
}  // namespace datadriven
}  // namespace sgpp
