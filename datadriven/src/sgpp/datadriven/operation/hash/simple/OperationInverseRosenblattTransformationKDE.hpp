// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// some defines for the following algorithm
#ifndef OPERATIONINVERSEROSENBLATTTRANSFORMATIONKDE_HPP
#define OPERATIONINVERSEROSENBLATTTRANSFORMATIONKDE_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/application/KernelDensityEstimator.hpp>

#include <sgpp/globaldef.hpp>

#include <random>

namespace sgpp {
namespace datadriven {

/**
 * Do inverse transformation in all dimensions
 */
class OperationInverseRosenblattTransformationKDE {
 public:
  OperationInverseRosenblattTransformationKDE(datadriven::KernelDensityEstimator& kde,
                                              double sigmaFactor = 6.0,
                                              double inversionEpsilon = 1e-10,
                                              std::uint64_t seed = std::mt19937_64::default_seed);
  virtual ~OperationInverseRosenblattTransformationKDE();

  /**
   * inverse Rosenblatt Transformation with mixed starting dimension
   *
   * @param pointsUniform data points to be transformed DataMatrix (rows: # of samples, columns: #
   * of dims)
   * @param pointsCdf Output DataMatrix (rows: # of samples, columns: # of dims)
   */
  virtual void doTransformation(base::DataMatrix& pointsUniform, base::DataMatrix& pointsCdf);

  virtual void doShuffledTransformation(base::DataMatrix& pointsUniform,
                                        base::DataMatrix& pointsCdf);

  /**
   * do the inverse Rosenblatt transformation for one data point for given samples
   *
   * @param y data point to be inverted
   * @param samples1d training samples in the dimension to be transformed
   * @param sigma bandwidth of the kernels in the current dimension
   * @param xlower lower bound for x-space
   * @param xupper upper bound for x-space
   * @param ylower lower bound for y-space
   * @param yupper upper bound for y-space
   * @param kern kernel evaluations
   * @return error of inversion
   */
  double doTransformation1D(double y, base::DataVector& samples1d, double sigma, double xlower,
                            double xupper, double ylower, double yupper, base::DataVector& kern);

  /// get the maximum error made during inversion
  double getMaxInversionError();

 private:
  datadriven::KernelDensityEstimator* kde;
  base::DataVector bandwidths;

  base::DataMatrix xlimits;
  base::DataMatrix ylimits;

  size_t ndim;
  size_t nsamples;

  /// maximum allowed inversion error
  double inversionEpsilon;

  /// shuffling devices
  std::mt19937_64 rng;

  /**
   * recalculates the search interval for bisection
   * @param sigmaFactor
   */
  void recalcLimits(double sigmaFactor);

  /**
   * Root finding using bisection algorithm for inverse CDF of KDE
   * @param y point where the CDF should be inverted
   * @param x root
   * @param xlower lower boundary of search interval
   * @param xupper upper boundary of search interval
   * @param samples1d samples in current dimension
   * @param sigma bandwidth for KDE
   * @param kern kernel evaluations for already processed dimensions
   * @param denom denominator for conditionalization of pdf
   * @param xacc accuracy
   * @param maxIterations maximum number of iterations
   */
  double bisection(double y, double& x, double& xlower, double& xupper, base::DataVector& samples1d,
                   double sigma, base::DataVector& kern, double denom, double xacc = 1e-8,
                   size_t maxIterations = 1000);

  /**
   * Root finding using newton's algorithm for inverse CDF of KDE
   * @param y point where the CDF should be inverted
   * @param x root
   * @param samples1d samples in current dimension
   * @param sigma bandwidth for KDE
   * @param kern kernel evaluations for already processed dimensions
   * @param denom denominator for conditionalization of pdf
   * @param xacc accuracy
   * @param maxIterations maximum number of iterations
   */
  double newton(double y, double& x, base::DataVector& samples1d, double sigma,
                base::DataVector& kern, double denom, double xacc = 1e-10,
                size_t maxIterations = 20);
};
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONINVERSEROSENBLATTTRANSFORMATIONKDE_HPP */
