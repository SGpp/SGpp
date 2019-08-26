// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/datadriven/application/DensityEstimator.hpp>
#include <sgpp/datadriven/application/KernelDensityEstimator.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationKDE.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationKDE.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalizeKDE.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/optimization/optimizer/unconstrained/NelderMead.hpp>
#include <sgpp/globaldef.hpp>

#include <map>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>

namespace sgpp {
namespace datadriven {

// -------------------- constructors and desctructors --------------------
KernelDensityEstimator::KernelDensityEstimator(KernelType kernelType,
                                               BandwidthOptimizationType bandwidthOptimizationType)
    : nsamples(0),
      ndim(0),
      bandwidths(0),
      norm(0),
      cond(0),
      sumCondInv(1.0),
      bandwidthOptimizationType(bandwidthOptimizationType) {
  initializeKernel(kernelType);
}

KernelDensityEstimator::KernelDensityEstimator(
    std::vector<std::shared_ptr<base::DataVector>>& samplesVec, KernelType kernelType,
    BandwidthOptimizationType bandwidthOptimizationType)
    : nsamples(0.0),
      ndim(samplesVec.size()),
      bandwidths(samplesVec.size()),
      norm(samplesVec.size()),
      cond(0.0),
      sumCondInv(0.0),
      bandwidthOptimizationType(bandwidthOptimizationType) {
  initializeKernel(kernelType);
  initialize(samplesVec);
}

KernelDensityEstimator::KernelDensityEstimator(base::DataMatrix& samples, KernelType kernelType,
                                               BandwidthOptimizationType bandwidthOptimizationType)
    : nsamples(samples.getNrows()),
      ndim(samples.getNcols()),
      bandwidths(samples.getNcols()),
      norm(samples.getNcols()),
      cond(samples.getNrows()),
      sumCondInv(0.0),
      bandwidthOptimizationType(bandwidthOptimizationType) {
  initializeKernel(kernelType);
  initialize(samples);
}

KernelDensityEstimator::KernelDensityEstimator(const KernelDensityEstimator& kde) {
  samplesVec = kde.samplesVec;
  nsamples = kde.nsamples;
  ndim = kde.ndim;
  bandwidths = base::DataVector(kde.bandwidths);
  norm = base::DataVector(kde.norm);
  cond = base::DataVector(kde.cond);
  sumCondInv = kde.sumCondInv;
  bandwidthOptimizationType = kde.bandwidthOptimizationType;

  initializeKernel(kde.kernel->getType());
}

KernelDensityEstimator::~KernelDensityEstimator() {}
// ----------------------------------------------------------------------

void KernelDensityEstimator::initializeKernel(KernelType kernelType) {
  switch (kernelType) {
    case KernelType::GAUSSIAN:
      kernel.reset(new GaussianKernel());
      break;
    case KernelType::EPANECHNIKOV:
      kernel.reset(new EpanechnikovKernel());
      break;
  }
}

void KernelDensityEstimator::initialize(base::DataMatrix& samples) {
  ndim = samples.getNcols();
  nsamples = samples.getNrows();

  samples.transpose();

  if (ndim > 0) {
    if (nsamples > 1) {
      // copy 1d samples to vector
      samplesVec.resize(ndim);

      for (size_t idim = 0; idim < ndim; idim++) {
        // copy
        samplesVec[idim] = std::make_shared<base::DataVector>(nsamples);
        samples.getRow(idim, *(samplesVec[idim]));
      }

      // initialize conditionalization factor
      cond.resize(nsamples);
      cond.setAll(1.0);
      sumCondInv = 1. / static_cast<double>(nsamples);

      // initialize normalization factors
      norm.resize(ndim);

      // init the bandwidths
      bandwidths.resize(ndim);
      computeAndSetOptKDEbdwth();
    } else {
      throw base::data_exception(
          "KernelDensityEstimator::KernelDensityEstimator: KDE needs at least two samples to "
          "estimate the bandwidth");
    }
  } else {
    throw base::data_exception(
        "KernelDensityEstimator::KernelDensityEstimator: KDE needs at least one dimensional data");
  }

  samples.transpose();
}

void KernelDensityEstimator::initialize(std::vector<std::shared_ptr<base::DataVector>>& samples) {
  ndim = samples.size();

  if (ndim > 0) {
    nsamples = samples[0]->getSize();

    if (nsamples > 0) {
      // copy 1d samples to vector
      samplesVec.resize(ndim);

      for (size_t idim = 0; idim < ndim; idim++) {
        samplesVec[idim] = std::make_shared<base::DataVector>(*(samples[idim]));  // copy
      }

      // initialize conditionalization factors
      cond.resize(nsamples);
      cond.setAll(1.0);
      sumCondInv = 1. / static_cast<double>(nsamples);

      // initialize normalization factors
      norm.resize(ndim);

      // init the bandwidths
      bandwidths.resize(ndim);
      computeAndSetOptKDEbdwth();
    } else {
      throw base::data_exception(
          "KernelDensityEstimator::KernelDensityEstimator : KDE needs at least two samples to "
          "estimate the bandwidth");
    }
  } else {
    throw base::data_exception(
        "KernelDensityEstimator::KernelDensityEstimator : KDE needs at least one dimensional data");
  }
}

size_t KernelDensityEstimator::getDim() { return ndim; }

size_t KernelDensityEstimator::getNsamples() { return nsamples; }

std::shared_ptr<base::DataMatrix> KernelDensityEstimator::getSamples() {
  std::shared_ptr<base::DataMatrix> ans = std::make_shared<base::DataMatrix>(ndim, nsamples);

  for (size_t idim = 0; idim < ndim; idim++) {
    ans->setRow(idim, *samplesVec[idim]);
  }

  ans->transpose();
  return ans;
}

std::shared_ptr<base::DataVector> KernelDensityEstimator::getSamples(size_t dim) {
  if (dim >= samplesVec.size()) {
    throw base::data_exception("KernelDensityEstimator::getSamples : dim out of range");
  }

  return samplesVec[dim];
}

void KernelDensityEstimator::getSample(size_t isample, base::DataVector& sample) {
  for (size_t idim = 0; idim < ndim; idim++) {
    sample[idim] = samplesVec[idim]->get(isample);
  }
}

void KernelDensityEstimator::getBandwidths(base::DataVector& sigma) {
  // copy
  sigma.resize(bandwidths.getSize());

  for (size_t i = 0; i < bandwidths.getSize(); i++) {
    sigma[i] = bandwidths[i];
  }
}

void KernelDensityEstimator::setBandwidths(const base::DataVector& sigma) {
  for (size_t i = 0; i < sigma.getSize(); i++) {
    bandwidths[i] = sigma[i];
    norm[i] = kernel->norm() / bandwidths[i];
  }
}

void KernelDensityEstimator::pdf(base::DataMatrix& data, base::DataVector& res) {
  // init variables
  base::DataVector x(ndim);

  // resize result vector
  res.resize(data.getNrows());
  res.setAll(0.0);

  // run over all data points
  for (size_t idata = 0; idata < data.getNrows(); idata++) {
    // copy samples
    for (size_t idim = 0; idim < ndim; idim++) {
      x[idim] = data.get(idata, idim);
    }

    res[idata] = pdf(x);
  }
}

double KernelDensityEstimator::pdf(base::DataVector& x) {
  // init variables
  double res = 0.0;

  // run over all data points
  for (size_t isample = 0; isample < nsamples; isample++) {
    res += evalKernel(x, isample);
  }

  return res * sumCondInv;
}

double KernelDensityEstimator::evalSubset(base::DataVector& x, std::vector<size_t> skipElements) {
  // init variables
  double res = 0.0;

  // sort the elements to be skipped
  std::sort(skipElements.begin(), skipElements.end());
  size_t isample = 0, j = 0;

  // just add those kernels which are not in the skipElements list
  while (isample < nsamples) {
    if (isample < skipElements[j]) {
      res += evalKernel(x, isample);
    } else {
      j++;
    }
    isample++;
  }

  return res / static_cast<double>(nsamples - skipElements.size());
}

double KernelDensityEstimator::evalKernel(base::DataVector& x, size_t i) {
  double res = 1.0;
  double y = 0.0;

  for (size_t idim = 0; idim < ndim; idim++) {
    // normalize x
    y = (x[idim] - samplesVec[idim]->get(i)) / bandwidths[idim];
    // evaluate kernel
    res *= norm[idim] * kernel->eval(y);
  }

  return cond[i] * res;
}

void KernelDensityEstimator::cov(base::DataMatrix& cov, base::DataMatrix* bounds) {
  if ((cov.getNrows() != ndim) || (cov.getNcols() != ndim)) {
    // covariance matrix has wrong size -> resize
    cov.resize(ndim, ndim);
  }

  // prepare covariance marix
  cov.setAll(0.0);

  // generate 1d densities and compute means and variances
  std::vector<double> means(ndim);
  std::vector<double> variances(ndim);

  std::unique_ptr<datadriven::OperationDensityMarginalizeKDE> opMarg(
      op_factory::createOperationDensityMarginalizeKDE(*this));
  KernelDensityEstimator kdeMarginalized(kernel->getType(), bandwidthOptimizationType);

  for (size_t idim = 0; idim < ndim; idim++) {
    opMarg->margToDimX(idim, kdeMarginalized);
    // store moments
    means[idim] = kdeMarginalized.mean();
    variances[idim] = kdeMarginalized.variance();
  }

  // helper variables
  std::vector<size_t> mdims(2);
  double covij = 0.0;

  KernelDensityEstimator kdeijdim(kernel->getType(), bandwidthOptimizationType);

  for (size_t idim = 0; idim < ndim; idim++) {
    // diagonal is equal to the variance of the marginalized densities
    cov.set(idim, idim, variances[idim]);

    for (size_t jdim = idim + 1; jdim < ndim; jdim++) {
      // marginalize the density
      mdims[0] = idim;
      mdims[1] = jdim;
      opMarg->margToDimXs(mdims, kdeijdim);
      // -----------------------------------------------------
      // compute the covariance of Cov(X_i, X_j)
      covij = kdeijdim.mean() - means[idim] * means[jdim];
      cov.set(idim, jdim, covij);
      cov.set(jdim, idim, covij);
      // -----------------------------------------------------
    }
  }
}

double KernelDensityEstimator::mean() {
  double res = 0, kernelMean = 1.;

  for (size_t isample = 0; isample < nsamples; isample++) {
    kernelMean = 1.;

    for (size_t idim = 0; idim < ndim; idim++) {
      kernelMean *= samplesVec[idim]->get(isample);
    }

    res += cond[isample] * kernelMean;
  }

  return res * sumCondInv;
}

double KernelDensityEstimator::variance() {
  double meansquared = 0, kernelVariance = 1., x = 0.0, sigma = 0.0;
  for (size_t isample = 0; isample < nsamples; isample++) {
    kernelVariance = 1.;
    for (size_t idim = 0; idim < ndim; idim++) {
      x = samplesVec[idim]->get(isample);
      kernelVariance *= sigma * sigma * kernel->variance() + x * x;
    }

    meansquared += cond[isample] * kernelVariance;
  }
  meansquared *= sumCondInv;

  double mu = mean();
  double var = meansquared - mu * mu;

  return var;
}

void KernelDensityEstimator::computeAndSetOptKDEbdwth() {
  base::DataVector sigma(ndim);

  switch (bandwidthOptimizationType) {
    case BandwidthOptimizationType::NONE:
      break;
    case BandwidthOptimizationType::SILVERMANSRULE:
      SilvermansRule::optimizeBandwidths(this, sigma);
      break;
    case BandwidthOptimizationType::SCOTTSRULE:
      ScottsRule::optimizeBandwidths(this, sigma);
      break;
    case BandwidthOptimizationType::MAXIMUMLIKELIHOOD:
      MaximumLikelihoodCrossValidation::optimizeBandwidths(this, sigma);
      break;
  }

  // set the current bandwidth
  setBandwidths(sigma);
}

// ------------------------- additional operations ---------------------------

void KernelDensityEstimator::getConditionalizationFactor(base::DataVector& pcond) {
  pcond.resize(nsamples);

  for (size_t isample = 0; isample < nsamples; isample++) {
    pcond[isample] = cond[isample];
  }
}

void KernelDensityEstimator::setConditionalizationFactor(base::DataVector& pcond) {
  double sumCond = 0.0;

  for (size_t isample = 0; isample < nsamples; isample++) {
    cond[isample] = pcond[isample];
    sumCond += cond[isample];
  }

  sumCondInv = 1. / sumCond;
}

void KernelDensityEstimator::updateConditionalizationFactors(base::DataVector& x,
                                                             std::vector<size_t>& dims,
                                                             base::DataVector& pcond) {
  // run over all samples and evaluate the kernels in each dimension
  // that should be conditionalized
  size_t idim = 0;
  double xi = 0.0;

  for (size_t i = 0; i < dims.size(); i++) {
    idim = dims[i];

    if (idim < ndim) {
      for (size_t isample = 0; isample < nsamples; isample++) {
        xi = (x[idim] - samplesVec[idim]->get(isample)) / bandwidths[idim];
        pcond[isample] *= norm[idim] * kernel->eval(xi);
      }
    } else {
      throw base::data_exception(
          "KernelDensityEstimator::updateConditionalizationFactors : can not conditionalize in non "
          "existing "
          "dimension");
    }
  }
}

Kernel& KernelDensityEstimator::getKernel() { return *kernel; }

KernelDensityEstimator* KernelDensityEstimator::margToDimX(size_t idim) {
  std::unique_ptr<datadriven::OperationDensityMarginalizeKDE> opMarg(
      op_factory::createOperationDensityMarginalizeKDE(*this));
  datadriven::KernelDensityEstimator* marginalizedKDE =
      new datadriven::KernelDensityEstimator(kernel->getType(), bandwidthOptimizationType);
  opMarg->margToDimX(idim, *marginalizedKDE);
  return marginalizedKDE;
}

KernelDensityEstimator* KernelDensityEstimator::marginalize(size_t idim) {
  std::unique_ptr<datadriven::OperationDensityMarginalizeKDE> opMarg(
      op_factory::createOperationDensityMarginalizeKDE(*this));
  datadriven::KernelDensityEstimator* marginalizedKDE =
      new datadriven::KernelDensityEstimator(kernel->getType(), bandwidthOptimizationType);
  opMarg->doMarginalize(idim, *marginalizedKDE);
  return marginalizedKDE;
}

// ----------------------------------------------------------------------------------
// kernels

Kernel::~Kernel() {}

GaussianKernel::~GaussianKernel() {}
double GaussianKernel::eval(double x) { return std::exp(-(x * x) / 2.); }
double GaussianKernel::cdf(double x) { return 0.5 + 0.5 * std::erf(x / M_SQRT2); }
double GaussianKernel::derivative(double x) { return x * eval(x); }
double GaussianKernel::norm() { return 1. / M_SQRT2PI; }
double GaussianKernel::variance() { return 1.0; }
KernelType GaussianKernel::getType() { return KernelType::GAUSSIAN; }

EpanechnikovKernel::~EpanechnikovKernel() {}

double EpanechnikovKernel::eval(double x) {
  if (x > -1 && x < 1.) {
    return 1. - x * x;
  } else {
    return 0.0;
  }
}

double EpanechnikovKernel::cdf(double x) {
  if (x < -1) {
    return 0.0;
  } else if (x < 1.) {
    return 0.75 * x * (1.0 - x * x / 3.) + 0.5;
  } else {
    return 1.0;
  }
}

double EpanechnikovKernel::derivative(double x) {
  if (x > -1 && x < 1.) {
    return 1.5 * x;
  } else {
    return 0.0;
  }
}

double EpanechnikovKernel::norm() { return 0.75; }
double EpanechnikovKernel::variance() { return 0.2; }

KernelType EpanechnikovKernel::getType() { return KernelType::EPANECHNIKOV; }

// ----------------------------------------------------------------------------------
// bandwidth optimizers

double RuleOfThumb::getSampleMean(base::DataVector& data) {
  double res = 0.;
  size_t n = data.getSize();

  for (size_t i = 0; i < n; i++) {
    res += data[i];
  }

  return res / static_cast<double>(n);
}

double RuleOfThumb::getSampleVariance(base::DataVector& data) {
  double mean = getSampleMean(data);
  double diff1 = 0.0;
  double diff2 = 0.0;

  size_t n = data.getSize();

  for (size_t i = 0; i < n; i++) {
    diff1 += (data[i] - mean) * (data[i] - mean);
    diff2 += (data[i] - mean);
  }

  return 1. / (static_cast<double>(n) - 1.) * (diff1 - 1. / static_cast<double>(n) * diff2 * diff2);
}

double RuleOfThumb::getSampleStd(base::DataVector& data) {
  return std::sqrt(getSampleVariance(data));
}
// ----------------------------------------------------------------------------------
void SilvermansRule::optimizeBandwidths(KernelDensityEstimator* kde, base::DataVector& bandwidths) {
  size_t numDims = kde->getDim();
  bandwidths.resize(numDims);

  base::DataVector flag(numDims);
  flag.setAll(1.);

  // get min and max in each direction
  double datamin = 0.0;
  double datamax = 0.0;
  std::shared_ptr<base::DataVector> samples1d;

  double stdd;

  for (size_t idim = 0; idim < numDims; idim++) {
    size_t numBorder = 0;
    samples1d = kde->getSamples(idim);
    size_t numSamples = samples1d->getSize();

    // search for maximum in current dimension
    datamin = samples1d->min();
    datamax = samples1d->max();

    double nearBorder = (datamax - datamin) / 20.;

    // count how many values are close to the border
    for (size_t isample = 0; isample < numSamples; isample++) {
      if (samples1d->get(isample) - datamin < nearBorder ||
          datamax - samples1d->get(isample) < nearBorder) {
        numBorder++;
      }
    }

    if (static_cast<double>(numBorder) > static_cast<double>(numSamples) / 20.) {
      flag[idim] = 0.5;
    }

    // compute the standard deviation
    stdd = RuleOfThumb::getSampleStd(*samples1d);

    // compute the bandwidth in dimension idim
    bandwidths[idim] =
        flag[idim] * std::pow(4. / (static_cast<double>(numDims) + 2),
                              1. / (static_cast<double>(numDims) + 4.)) *
        stdd * std::pow(static_cast<double>(numSamples), -1. / (static_cast<double>(numDims) + 4.));
  }
}
// ----------------------------------------------------------------------------------
void ScottsRule::optimizeBandwidths(KernelDensityEstimator* kde, base::DataVector& bandwidths) {
  size_t numDims = kde->getDim();
  bandwidths.resize(numDims);

  std::shared_ptr<base::DataVector> samples1d;
  double stdd;

  for (size_t idim = 0; idim < numDims; idim++) {
    samples1d = kde->getSamples(idim);
    size_t numSamples = samples1d->getSize();

    // compute the standard deviation
    stdd = RuleOfThumb::getSampleStd(*samples1d);

    // compute the bandwidth in dimension idim
    bandwidths[idim] =
        stdd * std::pow(static_cast<double>(numSamples), -1. / (static_cast<double>(numDims) + 4.));
  }
}
// ----------------------------------------------------------------------------------

KDEMaximumLikelihoodCrossValidation::KDEMaximumLikelihoodCrossValidation(
    KernelDensityEstimator& kde, size_t kfold, std::uint64_t seedValue)
    : sgpp::base::ScalarFunction(kde.getDim()), kde(kde), strain(kfold), stest(kfold) {
  // split the data set
  auto samples = kde.getSamples();
  size_t numSamples = samples->getNrows();
  size_t numDims = samples->getNcols();

  // choose randomly indices to get the training data set
  std::mt19937_64 gen(seedValue);
  std::uniform_int_distribution<size_t> dist(0, numSamples);

  base::DataVector p(numDims);
  base::DataVector tmp(numDims);

  std::vector<size_t> s(kfold);        // size of partition
  std::vector<size_t> ind(kfold + 1);  // index of partition

  // set size of partitions
  ind[0] = 0;
  for (size_t i = 0; i < kfold - 1; i++) {
    s[i] = numSamples / kfold;
    ind[i + 1] = ind[i] + s[i];
  }
  ind[kfold] = numSamples;
  s[kfold - 1] = numSamples - (kfold - 1) * (numSamples / kfold);

  // fill data
  for (size_t i = 0; i < kfold; i++) {
    // allocate memory
    strain[i] = std::make_shared<base::DataMatrix>(numSamples - s[i], numDims);
    stest[i] = std::make_shared<base::DataMatrix>(s[i], numDims);

    size_t localTest = 0;
    size_t localTrain = 0;

    for (size_t j = 0; j < numSamples; j++) {
      samples->getRow(j, p);

      if (ind[i] <= j && j < ind[i + 1]) {
        stest[i]->setRow(localTest, p);
        localTest++;
      } else {
        strain[i]->setRow(localTrain, p);
        localTrain++;
      }
    }
  }
}

void MaximumLikelihoodCrossValidation::optimizeBandwidths(KernelDensityEstimator* kde,
                                                          base::DataVector& bandwidths) {
  KDEMaximumLikelihoodCrossValidation f(*kde);
  sgpp::optimization::optimizer::NelderMead nelderMead(f, 200);
  nelderMead.optimize();
  base::DataVector xOptNM = nelderMead.getOptimalPoint();

  // copy result to bandwidhts
  for (size_t i = 0; i < xOptNM.getSize(); i++) {
    bandwidths[i] = xOptNM[i];
  }
}

double KDEMaximumLikelihoodCrossValidation::eval(const base::DataVector& x) {
  double result = 0.0;
  // do the k-fold cross validation
  for (size_t k = 0; k < strain.size(); k++) {
    // load the data set
    auto trainSamples = strain[k];
    auto testSamples = stest[k];

    KernelDensityEstimator localKDE(*trainSamples, kde.getKernel().getType(),
                                    BandwidthOptimizationType::NONE);
    localKDE.setBandwidths(x);

    // compute the cross entropy
    result += localKDE.crossEntropy(*testSamples);
  }

  return result / static_cast<double>(strain.size());

  //  // MLCV
  //  double result = 0.0;
  //
  //  // do the k-fold cross validation
  //  base::DataVector sample(kde.getDim());
  //  std::vector<size_t> skipElements(1);
  //  for (size_t k = 0; k < strain.size(); k++) {
  //    // load the data set
  //    auto trainSamples = strain[k];
  //    auto testSamples = stest[k];
  //
  //    KernelDensityEstimator localKDE(*trainSamples, kde.getKernel().getType(),
  //                                    BandwidthOptimizationType::NONE);
  //    localKDE.setBandwidths(x);
  //
  //    // compute the cross entropy
  //    double mlcv = 0.0;
  //    double value = 0.0;
  //    std::uint32_t count = 0;
  //    for (size_t i = 0; i < trainSamples->getNrows(); i++) {
  //      skipElements[0] = i;
  //      trainSamples->getRow(i, sample);
  //      value = kde.evalSubset(sample, skipElements);
  //      if (value > 1e-10) {
  //        mlcv -= std::log2(value);
  //        count += 1;
  //      }
  //    }
  //
  //    // update the result
  //    result += mlcv / static_cast<double>(count);
  //  }
  //
  //  return result / static_cast<double>(strain.size());
}

}  // namespace datadriven
}  // namespace sgpp
