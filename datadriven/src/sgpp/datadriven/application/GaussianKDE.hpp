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

namespace SGPP {
namespace datadriven {

#ifndef M_SQRT2PI
#define M_SQRT2PI       2.506628274631000241612355239340    /* sqrt(2*pi) */
#endif

class GaussianKDE: public DensityEstimator {
 public:
  GaussianKDE();
  GaussianKDE(std::vector<base::DataVector*>& samplesVec);
  GaussianKDE(base::DataMatrix& samples);
  virtual ~GaussianKDE();

  virtual void initialize(base::DataMatrix& samples);
  virtual void initialize(std::vector<base::DataVector*>& samplesVec);

  float_t mean();
  float_t variance();
  float_t std_deviation();

  void cov(base::DataMatrix& cov);

  float_t pdf(base::DataVector& x);
  void pdf(base::DataMatrix& points, base::DataVector& res);

  /// getter and setter functions
  void getConditionalizationFactor(base::DataVector& pcond);
  void setConditionalizationFactor(base::DataVector& pcond);
  void updateConditionalizationFactors(base::DataVector& x,
                                       std::vector<size_t>& dims, base::DataVector& pcond);

  void getBandwidths(base::DataVector& sigma);

  virtual base::DataMatrix* getSamples();
  virtual base::DataVector* getSamples(size_t dim);

  size_t getDim();
  size_t getNsamples();

 private:
  /// samples
  std::vector<base::DataVector*> samplesVec;

  size_t nsamples;
  size_t ndim;

  /// standard deviations for the kernels in 1d
  base::DataVector bandwidths;
  /// normalization factor for 1d kernels
  base::DataVector norm;
  /// conditionalization factors
  base::DataVector cond;
  float_t sumCond;

  void computeOptKDEbdwth();
  void computeNormalizationFactors();

  float_t getSampleMean(base::DataVector& data);
  float_t getSampleVariance(base::DataVector& data);
  float_t getSampleStd(base::DataVector& data);
};

}  // namespace datadriven
}  // namespace sg

#endif /* GAUSSIANGAUSSIANKDE_HPP_ */
