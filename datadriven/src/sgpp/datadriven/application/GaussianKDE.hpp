// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GAUSSIANGAUSSIANKDE_HPP_
#define GAUSSIANGAUSSIANKDE_HPP_

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/application/DensityEstimator.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

#ifndef M_SQRT2PI
#define M_SQRT2PI 2.506628274631000241612355239340 /* sqrt(2*pi) */
#endif

class GaussianKDE : public DensityEstimator {
 public:
  GaussianKDE();
  explicit GaussianKDE(std::vector<std::shared_ptr<base::DataVector>>& samplesVec);
  explicit GaussianKDE(base::DataMatrix& samples);
  virtual ~GaussianKDE();

  virtual void initialize(base::DataMatrix& samples);
  virtual void initialize(std::vector<std::shared_ptr<base::DataVector>>& samplesVec);

  double mean();
  double variance();
  double std_deviation();

  void cov(base::DataMatrix& cov);

  double pdf(base::DataVector& x);
  void pdf(base::DataMatrix& points, base::DataVector& res);

  /// getter and setter functions
  void getConditionalizationFactor(base::DataVector& pcond);
  void setConditionalizationFactor(base::DataVector& pcond);
  void updateConditionalizationFactors(base::DataVector& x, std::vector<size_t>& dims,
                                       base::DataVector& pcond);

  void getBandwidths(base::DataVector& sigma);

  virtual std::shared_ptr<base::DataMatrix> getSamples();
  virtual std::shared_ptr<base::DataVector> getSamples(size_t dim);

  size_t getDim();
  size_t getNsamples();

 private:
  /// samples
  std::vector<std::shared_ptr<base::DataVector>> samplesVec;

  size_t nsamples;
  size_t ndim;

  /// standard deviations for the kernels in 1d
  base::DataVector bandwidths;
  /// normalization factor for 1d kernels
  base::DataVector norm;
  /// conditionalization factors
  base::DataVector cond;
  double sumCond;

  void computeOptKDEbdwth();
  void computeNormalizationFactors();

  double getSampleMean(base::DataVector& data);
  double getSampleVariance(base::DataVector& data);
  double getSampleStd(base::DataVector& data);
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* GAUSSIANGAUSSIANKDE_HPP_ */
