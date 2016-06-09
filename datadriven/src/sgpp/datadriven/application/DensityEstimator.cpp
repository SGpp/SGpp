// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/DensityEstimator.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <algorithm>

namespace sgpp {
namespace datadriven {

// -------------------- constructors and desctructors --------------------
DensityEstimator::DensityEstimator() {}

DensityEstimator::~DensityEstimator() {}
// ----------------------------------------------------------------------

double DensityEstimator::std_deviation() { return std::sqrt(variance()); }

void DensityEstimator::corrcoef(base::DataMatrix& corr) {
  // get covariance matrix and ...
  cov(corr);

  // ... normalize it
  double corrij = 0.0;
  double sigmai = 0.0, sigmaj = 0.0;
  size_t ndim = corr.getNcols();

  for (size_t idim = 0; idim < ndim; idim++) {
    // normalize row idim but the diagonal element
    sigmai = sqrt(corr.get(idim, idim));

    for (size_t jdim = idim + 1; jdim < ndim; jdim++) {
      sigmaj = sqrt(corr.get(jdim, jdim));
      corrij = corr.get(idim, jdim) / (sigmai * sigmaj);
      corr.set(idim, jdim, corrij);
      corr.set(jdim, idim, corrij);
    }

    // set the diagonal element
    corr.set(idim, idim, 1.0);
  }
}

double DensityEstimator::crossEntropy(sgpp::base::DataMatrix& samples) {
  size_t numSamples = samples.getNrows();
  size_t numDims = samples.getNcols();

  if (numSamples > 0) {
    base::DataVector sample(numDims);
    double sum = 0.0;
    for (size_t i = 0; i < numSamples; i++) {
      samples.getRow(i, sample);
      sum += std::log2(std::max(1e-10, pdf(sample)));
    }

    return -1.0 * sum / static_cast<double>(numSamples);
  } else {
    throw sgpp::base::algorithm_exception(
        "DensityEstimator::crossEntropy - size of test samples is zero");
  }
}

}  // namespace datadriven
}  // namespace sgpp
