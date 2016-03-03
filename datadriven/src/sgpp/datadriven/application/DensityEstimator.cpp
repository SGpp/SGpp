// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/DensityEstimator.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace SGPP {
namespace datadriven {

// -------------------- constructors and desctructors --------------------
DensityEstimator::DensityEstimator() {}

DensityEstimator::~DensityEstimator() {}
// ----------------------------------------------------------------------

float_t DensityEstimator::std_deviation() { return std::sqrt(variance()); }

void DensityEstimator::corrcoef(base::DataMatrix& corr) {
  // get covariance matrix and ...
  cov(corr);

  // ... normalize it
  float_t corrij = 0.0;
  float_t sigmai = 0.0, sigmaj = 0.0;
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

}  // namespace datadriven
}  // namespace SGPP
