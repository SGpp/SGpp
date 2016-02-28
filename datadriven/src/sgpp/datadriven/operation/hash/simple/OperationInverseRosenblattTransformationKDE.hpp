/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de
// some defines for the following algorithm
#ifndef OPERATIONINVERSEROSENBLATTTRANSFORMATIONKDE_HPP
#define OPERATIONINVERSEROSENBLATTTRANSFORMATIONKDE_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/application/GaussianKDE.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

/**
 * Do inverse transformation in all dimensions
 */
class OperationInverseRosenblattTransformationKDE {
 public:
  OperationInverseRosenblattTransformationKDE(datadriven::GaussianKDE& kde,
                                              float_t sigmaFactor = 6.0,
                                              float_t inversionEpsilon = 1e-10);
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
  float_t doTransformation1D(float_t y, base::DataVector& samples1d, float_t sigma, float_t xlower,
                             float_t xupper, float_t ylower, float_t yupper,
                             base::DataVector& kern);

  /// get the maximum error made during inversion
  float_t getMaxInversionError();

 private:
  datadriven::GaussianKDE* kde;
  base::DataVector bandwidths;

  base::DataMatrix xlimits;
  base::DataMatrix ylimits;

  size_t ndim;
  size_t nsamples;

  /// maximum allowed inversion error
  float_t inversionEpsilon;

  /**
   * recalculates the search interval for bisection
   * @param sigmaFactor
   */
  void recalcLimits(float_t sigmaFactor);

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
  float_t bisection(float_t y, float_t& x, float_t& xlower, float_t& xupper,
                    base::DataVector& samples1d, float_t sigma, base::DataVector& kern,
                    float_t denom, float_t xacc = 1e-8, size_t maxIterations = 1000);

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
  float_t newton(float_t y, float_t& x, base::DataVector& samples1d, float_t sigma,
                 base::DataVector& kern, float_t denom, float_t xacc = 1e-10,
                 size_t maxIterations = 20);

  /**
   * Root finding using halley's algorithm for inverse CDF of KDE
   * @param y point where the CDF should be inverted
   * @param x root
   * @param samples1d samples in current dimension
   * @param sigma bandwidth for KDE
   * @param kern kernel evaluations for already processed dimensions
   * @param denom denominator for conditionalization of pdf
   * @param xacc accuracy
   * @param maxIterations maximum number of iterations
   */
  float_t halley(float_t y, float_t& x, base::DataVector& samples1d, float_t sigma,
                 base::DataVector& kern, float_t denom, float_t xacc = 1e-10,
                 size_t maxIterations = 20);
};
}  // namespace datadriven
}  // namespace SGPP
#endif /* OPERATIONINVERSEROSENBLATTTRANSFORMATIONKDE_HPP */
