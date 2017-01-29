// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/DensityEstimator.hpp>
#include <sgpp/datadriven/application/GaussianKDE.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationKDE.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationKDE.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalizeKDE.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/globaldef.hpp>

#include <map>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <random>
#include <vector>

namespace sgpp {
namespace datadriven {

// -------------------- constructors and desctructors --------------------
GaussianKDE::GaussianKDE() : nsamples(0), ndim(0), bandwidths(0), norm(0), cond(0), sumCond(1.0) {}

GaussianKDE::GaussianKDE(std::vector<std::shared_ptr<base::DataVector>>& samplesVec)
    : nsamples(0.0),
      ndim(samplesVec.size()),
      bandwidths(samplesVec.size()),
      norm(samplesVec.size()),
      cond(0.0),
      sumCond(0.0) {
  initialize(samplesVec);
}

GaussianKDE::GaussianKDE(base::DataMatrix& samples)
    : nsamples(samples.getNrows()),
      ndim(samples.getNcols()),
      bandwidths(samples.getNcols()),
      norm(samples.getNcols()),
      cond(samples.getNrows()),
      sumCond(0.0) {
  initialize(samples);
}

GaussianKDE::GaussianKDE(const GaussianKDE& kde) {
  samplesVec = kde.samplesVec;
  nsamples = kde.nsamples;
  ndim = kde.ndim;
  bandwidths = base::DataVector(kde.bandwidths);
  norm = base::DataVector(kde.norm);
  cond = base::DataVector(kde.cond);
  sumCond = kde.sumCond;
}

GaussianKDE::~GaussianKDE() {}
// ----------------------------------------------------------------------

void GaussianKDE::initialize(base::DataMatrix& samples) {
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

      // init the bandwidths
      bandwidths.resize(ndim);
      computeOptKDEbdwth();

      // initialize normalization factors
      norm.resize(ndim);

      for (size_t d = 0; d < ndim; d++) {
        norm[d] = 1. / (bandwidths[d] * M_SQRT2PI);
      }

      // initialize conditionalization factor
      cond.resize(nsamples);
      cond.setAll(1.0);
      sumCond = static_cast<double>(nsamples);
    } else {
      throw base::data_exception(
          "GaussianKDE::GaussianKDE: KDE needs at least two samples to estimate the bandwidth");
    }
  } else {
    throw base::data_exception("GaussianKDE::GaussianKDE: KDE needs at least one dimensional data");
  }

  samples.transpose();
}

void GaussianKDE::initialize(std::vector<std::shared_ptr<base::DataVector>>& samples) {
  ndim = samples.size();

  if (ndim > 0) {
    nsamples = samples[0]->getSize();

    if (nsamples > 0) {
      // copy 1d samples to vector
      samplesVec.resize(ndim);

      for (size_t idim = 0; idim < ndim; idim++) {
        samplesVec[idim] = std::make_shared<base::DataVector>(*(samples[idim]));  // copy
      }

      // init the bandwidths
      bandwidths.resize(ndim);
      computeOptKDEbdwth();

      // initialize normalization factors
      norm.resize(ndim);

      for (size_t d = 0; d < ndim; d++) {
        norm[d] = 1. / (bandwidths[d] * M_SQRT2PI);
      }

      // initialize conditionalization factors
      cond.resize(nsamples);
      cond.setAll(1.0);
      sumCond = static_cast<double>(nsamples);
    } else {
      throw base::data_exception(
          "GaussianKDE::GaussianKDE : KDE needs at least two samples to estimate the bandwidth");
    }
  } else {
    throw base::data_exception(
        "GaussianKDE::GaussianKDE : KDE needs at least one dimensional data");
  }
}

size_t GaussianKDE::getDim() { return ndim; }

size_t GaussianKDE::getNsamples() { return nsamples; }

std::shared_ptr<base::DataMatrix> GaussianKDE::getSamples() {
  std::shared_ptr<base::DataMatrix> ans = std::make_shared<base::DataMatrix>(ndim, nsamples);

  for (size_t idim = 0; idim < ndim; idim++) {
    ans->setRow(idim, *samplesVec[idim]);
  }

  ans->transpose();
  return ans;
}

std::shared_ptr<base::DataVector> GaussianKDE::getSamples(size_t dim) {
  if (dim >= samplesVec.size()) {
    throw base::data_exception("GaussianKDE::getSamples : dim out of range");
  }

  return samplesVec[dim];
}

void GaussianKDE::getBandwidths(base::DataVector& sigma) {
  // copy
  sigma.resize(bandwidths.getSize());

  for (size_t i = 0; i < bandwidths.getSize(); i++) {
    sigma[i] = bandwidths[i];
  }
}

void GaussianKDE::pdf(base::DataMatrix& data, base::DataVector& res) {
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

double GaussianKDE::pdf(base::DataVector& x) {
  // init variables
  double res = 0.0;
  double kern = 0, y = 0.0;

  // run over all data points
  for (size_t isample = 0; isample < nsamples; isample++) {
    kern = 1.;

    for (size_t idim = 0; idim < ndim; idim++) {
      // normalize x
      y = (x[idim] - samplesVec[idim]->get(isample)) / bandwidths[idim];
      // evaluate kernel
      kern *= norm[idim] * std::exp(-(y * y) / 2.);
    }

    res += cond[isample] * kern;
  }

  return res / sumCond;
}

void GaussianKDE::cov(base::DataMatrix& cov) {
  if ((cov.getNrows() != ndim) || (cov.getNcols() != ndim)) {
    // throw error -> covariance matrix has wrong size
    throw base::data_exception("GaussianKDE::cov : covariance matrix has the wrong size");
  }

  // prepare covariance marix
  cov.setAll(0.0);

  // generate 1d densities and compute means and variances
  std::vector<double> means(ndim);
  std::vector<double> variances(ndim);

  std::unique_ptr<datadriven::OperationDensityMarginalizeKDE> opMarg(
      op_factory::createOperationDensityMarginalizeKDE(*this));
  GaussianKDE kdeMarginalized;

  for (size_t idim = 0; idim < ndim; idim++) {
    opMarg->margToDimX(idim, kdeMarginalized);
    // store moments
    means[idim] = kdeMarginalized.mean();
    variances[idim] = kdeMarginalized.variance();
  }

  // helper variables
  std::vector<size_t> mdims(2);
  double covij = 0.0;

  GaussianKDE kdeijdim;

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

double GaussianKDE::mean() {
  double res = 0, kernelMean = 1.;

  for (size_t isample = 0; isample < nsamples; isample++) {
    kernelMean = 1.;

    for (size_t idim = 0; idim < ndim; idim++) {
      kernelMean *= samplesVec[idim]->get(isample);
    }

    res += kernelMean;
  }

  return res / static_cast<double>(nsamples);
}

double GaussianKDE::variance() {
  double meansquared = 0, kernelVariance = 1., x = 0.0, sigma = 0.0;

  for (size_t isample = 0; isample < nsamples; isample++) {
    kernelVariance = 1.;

    for (size_t idim = 0; idim < ndim; idim++) {
      x = samplesVec[idim]->get(isample);
      sigma = bandwidths[idim];
      kernelVariance *= sigma * sigma + x * x;
    }

    meansquared += kernelVariance;
  }

  meansquared /= static_cast<double>(nsamples);

  double mu = mean();
  double var = meansquared - mu * mu;

  return var;
}

void GaussianKDE::computeOptKDEbdwth() {
  if (ndim != bandwidths.getSize()) {
    throw base::data_exception("GaussianKDE::computeOptKDEbdwth : KDEBdwth dimension error");
  }

  base::DataVector flag(ndim);
  flag.setAll(1.);

  // get min and max in each direction
  double datamin = 0.0;
  double datamax = 0.0;
  std::shared_ptr<base::DataVector> samples1d;

  double stdd;

  for (size_t idim = 0; idim < ndim; idim++) {
    size_t numBorder = 0;
    samples1d = samplesVec[idim];
    // search for maximum in current dimension
    datamin = samples1d->min();
    datamax = samples1d->max();

    double nearBorder = (datamax - datamin) / 20.;

    // count how many values are close to the border
    for (size_t isample = 0; isample < nsamples; isample++) {
      if (samples1d->get(isample) - datamin < nearBorder ||
          datamax - samples1d->get(isample) < nearBorder) {
        numBorder++;
      }
    }

    if (static_cast<double>(numBorder) > static_cast<double>(nsamples) / 20.) {
      flag[idim] = 0.5;
    }

    // compute the standard deviation
    stdd = getSampleStd(*samples1d);

    // compute the bandwidth in dimension idim
    bandwidths[idim] =
        flag[idim] *
        std::pow(4. / (static_cast<double>(ndim) + 2), 1. / (static_cast<double>(ndim) + 4.)) *
        stdd * std::pow(static_cast<double>(nsamples), -1. / (static_cast<double>(ndim) + 4.));
  }

  return;
}

double GaussianKDE::getSampleMean(base::DataVector& data) {
  double res = 0.;
  size_t n = data.getSize();

  for (size_t i = 0; i < n; i++) {
    res += data[i];
  }

  return res / static_cast<double>(n);
}

double GaussianKDE::getSampleVariance(base::DataVector& data) {
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

double GaussianKDE::getSampleStd(base::DataVector& data) {
  return std::sqrt(getSampleVariance(data));
}

// ------------------------- additional operations ---------------------------

void GaussianKDE::getConditionalizationFactor(base::DataVector& pcond) {
  pcond.resize(nsamples);

  for (size_t isample = 0; isample < nsamples; isample++) {
    pcond[isample] = cond[isample];
  }
}

void GaussianKDE::setConditionalizationFactor(base::DataVector& pcond) {
  sumCond = 0.0;

  for (size_t isample = 0; isample < nsamples; isample++) {
    cond[isample] = pcond[isample];
    sumCond += cond[isample];
  }
}

void GaussianKDE::updateConditionalizationFactors(base::DataVector& x, std::vector<size_t>& dims,
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
        pcond[isample] *= norm[idim] * std::exp(-(xi * xi) / 2.);
      }
    } else {
      throw base::data_exception(
          "GaussianKDE::updateConditionalizationFactors : can not conditionalize in non existing "
          "dimension");
    }
  }
}

}  // namespace datadriven
}  // namespace sgpp
