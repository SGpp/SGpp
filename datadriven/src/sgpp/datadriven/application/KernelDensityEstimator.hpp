// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/application/DensityEstimator.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

// --------------------------------------------------------------------------------
enum class KernelType { GAUSSIAN, EPANECHNIKOV };

enum class BandwidthOptimizationType { NONE, RULEOFTHUMB, MAXIMUMLIKELIHOOD };

#ifndef M_SQRT2PI
#define M_SQRT2PI 2.506628274631000241612355239340 /* sqrt(2*pi) */
#endif

class Kernel {
 public:
  virtual ~Kernel();
  virtual double eval(double x) = 0;
  virtual double cdf(double x) = 0;
  virtual double derivative(double x) = 0;

  virtual double norm() = 0;
  virtual KernelType getType() = 0;
};

class GaussianKernel : public Kernel {
 public:
  virtual ~GaussianKernel();

  double eval(double x) override;
  double cdf(double x) override;
  double derivative(double x) override;
  double norm() override;
  KernelType getType() override;
};

class EpanechnikovKernel : public Kernel {
 public:
  virtual ~EpanechnikovKernel();

  double eval(double x) override;
  double cdf(double x) override;
  double derivative(double x) override;
  double norm() override;
  KernelType getType() override;
};

// --------------------------------------------------------------------------------

class KernelDensityEstimator : public DensityEstimator {
 public:
  explicit KernelDensityEstimator(
      KernelType kernelType = KernelType::GAUSSIAN,
      BandwidthOptimizationType bandwidthOptimizationType = BandwidthOptimizationType::RULEOFTHUMB);
  explicit KernelDensityEstimator(
      std::vector<std::shared_ptr<base::DataVector>>& samplesVec,
      KernelType kernelType = KernelType::GAUSSIAN,
      BandwidthOptimizationType bandwidthOptimizationType = BandwidthOptimizationType::RULEOFTHUMB);
  explicit KernelDensityEstimator(
      base::DataMatrix& samples, KernelType kernelType = KernelType::GAUSSIAN,
      BandwidthOptimizationType bandwidthOptimizationType = BandwidthOptimizationType::RULEOFTHUMB);
  KernelDensityEstimator(const KernelDensityEstimator& kde);

  virtual ~KernelDensityEstimator();

  void initialize(base::DataMatrix& samples) override;
  void initialize(std::vector<std::shared_ptr<base::DataVector>>& samplesVec);
  void initializeKernel(KernelType kernelType);

  Kernel& getKernel();
  double mean();
  double variance();

  void cov(base::DataMatrix& cov);

  double pdf(base::DataVector& x);
  void pdf(base::DataMatrix& points, base::DataVector& res);

  double evalSubset(base::DataVector& x, std::vector<size_t> skipElements);

  /// getter and setter functions
  void getConditionalizationFactor(base::DataVector& pcond);
  void setConditionalizationFactor(base::DataVector& pcond);
  void updateConditionalizationFactors(base::DataVector& x, std::vector<size_t>& dims,
                                       base::DataVector& pcond);

  void getBandwidths(base::DataVector& sigma);
  void setBandwidths(const base::DataVector& sigma);

  std::shared_ptr<base::DataMatrix> getSamples() override;
  std::shared_ptr<base::DataVector> getSamples(size_t dim) override;
  void getSample(size_t isample, base::DataVector& sample);

  size_t getDim();
  size_t getNsamples();

 private:
  double evalKernel(base::DataVector& x, size_t i);

  /// samples
  std::vector<std::shared_ptr<base::DataVector>> samplesVec;

  /// kernel
  std::unique_ptr<Kernel> kernel;

  size_t nsamples;
  size_t ndim;

  /// standard deviations for the kernels in 1d
  base::DataVector bandwidths;
  /// normalization factor for 1d kernels
  base::DataVector norm;
  /// conditionalization factors
  base::DataVector cond;
  double sumCond;

  /// bandwith optimization type
  BandwidthOptimizationType bandwidthOptimizationType;

  void computeAndSetOptKDEbdwth();
  void computeNormalizationFactors();
};

// --------------------------------------------------------------------------------
class KDEMaximumLikelihoodCrossValidation : public sgpp::optimization::ScalarFunction {
 public:
  /**
   * Constructor.
   */
  explicit KDEMaximumLikelihoodCrossValidation(
      KernelDensityEstimator& kde, size_t kfold = 10,
      std::uint64_t seedValue = std::mt19937_64::default_seed);

  double eval(const sgpp::base::DataVector& x);

  /**
   * @param[out] clone pointer to cloned object
   */
  virtual void clone(std::unique_ptr<sgpp::optimization::ScalarFunction>& clone) const {
    clone = std::unique_ptr<sgpp::optimization::ScalarFunction>(
        new KDEMaximumLikelihoodCrossValidation(*this));
  }

 private:
  KernelDensityEstimator& kde;
  std::vector<std::shared_ptr<base::DataMatrix>> strain;
  std::vector<std::shared_ptr<base::DataMatrix>> stest;
};

// --------------------------------------------------------------------------------

class RuleOfThumb {
 public:
  virtual ~RuleOfThumb();
  static void optimizeBandwidths(KernelDensityEstimator* kde, base::DataVector& bandwidths);

 private:
  static double getSampleMean(base::DataVector& data);
  static double getSampleVariance(base::DataVector& data);
  static double getSampleStd(base::DataVector& data);
};

class MaximumLikelihoodCrossValidation {
 public:
  virtual ~MaximumLikelihoodCrossValidation();
  static void optimizeBandwidths(KernelDensityEstimator* kde, base::DataVector& bandwidths);
};

} /* namespace datadriven */
} /* namespace sgpp */
