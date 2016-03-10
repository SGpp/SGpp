// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONROSENBLATTTRANSFORMATIONKDE_HPP
#define OPERATIONROSENBLATTTRANSFORMATIONKDE_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/application/GaussianKDE.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Do transformation in all dimensions
 */
class OperationRosenblattTransformationKDE {
 public:
  explicit OperationRosenblattTransformationKDE(GaussianKDE& kde);
  virtual ~OperationRosenblattTransformationKDE();

  /**
   * Rosenblatt Transformation with mixed starting dimension
   *
   * @param pointsCdf Output base::DataMatrix (rows: # of samples, columns: # of dims)
   * @param pointsUniform data points to be transformed base::DataMatrix (rows: # of samples,
   * columns: # of dims)
   */
  virtual void doTransformation(base::DataMatrix& pointsCdf, base::DataMatrix& pointsUniform);

  virtual void doShuffledTransformation(base::DataMatrix& pointsCdf,
                                        base::DataMatrix& pointsUniform);

  /**
   * Rosenblatt transformation for one data point with given samples and
   * and kernel evaluations, see doTransformation for details.
   *
   * @param x data point
   * @param samples1d training samples in the dimension to be transformed
   * @param sigma bandwidth of the kernels in the current dimension
   * @param kern kernel evaluations
   */
  double doTransformation1D(double x, base::DataVector& samples1d, double sigma,
                             base::DataVector& kern);

 private:
  datadriven::GaussianKDE* kde;
  base::DataVector bandwidths;

  size_t ndim;
  size_t nsamples;
};

}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONINVERSEROSENBLATTTRANSFORMATIONKDE_HPP */
