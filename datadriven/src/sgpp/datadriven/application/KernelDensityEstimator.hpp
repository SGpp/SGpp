// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/application/DensityEstimator.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

enum class KernelType { GAUSSIAN, EPANECHNIKOV };

// enum class BandwithOptimizationType {
//  RULEOFTHUMB,
//  MLCV  // maximum likelihood cross-validation
//};

#ifndef M_SQRT2PI
#define M_SQRT2PI 2.506628274631000241612355239340 /* sqrt(2*pi) */
#endif

class Kernel {
 public:
  virtual ~Kernel() {}
  virtual double eval(double x) = 0;
  virtual double cdf(double x) = 0;
  virtual double derivative(double x) = 0;

  virtual double norm() = 0;
  virtual KernelType getType() = 0;
};

class GaussianKernel : public Kernel {
 public:
  virtual ~GaussianKernel() {}
  double eval(double x) override { return std::exp(-(x * x) / 2.); }
  double cdf(double x) override { return 0.5 + 0.5 * std::erf(x / M_SQRT2); }
  double derivative(double x) override { return x * eval(x); }  // / sigma

  double norm() override { return 1. / M_SQRT2PI; }

  KernelType getType() { return KernelType::GAUSSIAN; }
};

class EpanechnikovKernel : public Kernel {
 public:
  virtual ~EpanechnikovKernel() {}

  double eval(double x) override {
    if (x > -1 && x < 1.) {
      return 1. - x * x;
    } else {
      return 0.0;
    }
  }

  double cdf(double x) override {
    if (x < -1) {
      return 0.0;
    } else if (x < 1.) {
      return 0.75 * x * (1.0 - x * x / 3.) + 0.5;
    } else {
      return 1.0;
    }
  }

  double derivative(double x) override {
    if (x > -1 && x < 1.) {
      return 1.5 * x;
    } else {
      return 0.0;
    }
  }

  double norm() override { return 0.75; }

  KernelType getType() { return KernelType::EPANECHNIKOV; }
};

class KernelDensityEstimator : public DensityEstimator {
 public:
  explicit KernelDensityEstimator(KernelType kernelType = KernelType::GAUSSIAN);
  explicit KernelDensityEstimator(std::vector<std::shared_ptr<base::DataVector>>& samplesVec,
                                  KernelType kernelType = KernelType::GAUSSIAN);
  explicit KernelDensityEstimator(base::DataMatrix& samples,
                                  KernelType kernelType = KernelType::GAUSSIAN);
  KernelDensityEstimator(const KernelDensityEstimator& kde);

  virtual ~KernelDensityEstimator();

  virtual void initialize(base::DataMatrix& samples);
  virtual void initialize(std::vector<std::shared_ptr<base::DataVector>>& samplesVec);
  virtual void initializeKernel(KernelType kernelType);

  Kernel& getKernel();
  double mean();
  double variance();

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

  /// kernel
  Kernel* kernel;

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

} /* namespace datadriven */
} /* namespace sgpp */
