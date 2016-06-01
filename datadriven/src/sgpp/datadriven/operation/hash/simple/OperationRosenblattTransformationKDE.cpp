// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationKDE.hpp>

#include <sgpp/globaldef.hpp>
#include <map>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

// #define DEBUG_ROSENBLATT

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven {

OperationRosenblattTransformationKDE::OperationRosenblattTransformationKDE(
    KernelDensityEstimator& density)
    : kde(&density),
      bandwidths(density.getDim()),
      ndim(density.getDim()),
      nsamples(density.getNsamples()) {
  // get the bandwidth from the density for optimized
  density.getBandwidths(bandwidths);
}

OperationRosenblattTransformationKDE::~OperationRosenblattTransformationKDE() {}

void OperationRosenblattTransformationKDE::doTransformation(DataMatrix& pointsCdf,
                                                            DataMatrix& pointsUniform) {
#pragma omp parallel
  {
#pragma omp for schedule(dynamic)
    for (size_t idata = 0; idata < pointsCdf.getNrows(); idata++) {
      // Work arrays
      DataVector unif(ndim);
      DataVector cdf(ndim);
      DataVector kern(nsamples);
      std::shared_ptr<base::DataVector> samples1d;

      kern.setAll(1.0);
      pointsCdf.getRow(idata, cdf);

      for (size_t idim = 0; idim < ndim; idim++) {
        // get samples in current dimension
        samples1d = kde->getSamples(idim);

        // transform the point in the current dimension
        unif[idim] = doTransformation1D(cdf[idim], *samples1d, bandwidths[idim], kern);

        // Update the kernel for the next dimension
        for (size_t isamples = 0; isamples < nsamples; isamples++) {
          double xi = (cdf[idim] - samples1d->get(isamples)) / bandwidths[idim];
          kern[isamples] *=
              kde->getKernel().eval(xi);  // std::exp(-(xi * xi) / 2.);  (bw*sqrt(2*PI)) cancels;
        }
      }

      // write them to the output
      pointsUniform.setRow(idata, unif);
    }
  }
  return;
}

void OperationRosenblattTransformationKDE::doShuffledTransformation(DataMatrix& pointsCdf,
                                                                    DataMatrix& pointsUniform) {
  // Work arrays
  DataVector unif(ndim);
  DataVector cdf(ndim);
  DataVector kern(nsamples);
  std::shared_ptr<base::DataVector> samples1d;

  double xi = 0;

  std::vector<size_t> dims(ndim);

  for (size_t i = 0; i < ndim; i++) {
    dims[i] = i;
  }

  for (size_t idata = 0; idata < pointsCdf.getNrows(); idata++) {
    kern.setAll(1.0);
    pointsCdf.getRow(idata, cdf);

    std::next_permutation(dims.begin(), dims.end());

    for (size_t i = 0; i < ndim; i++) {
      size_t idim = dims[i];
      // get samples in current dimension
      samples1d = kde->getSamples(idim);

      // transform the point in the current dimension
      unif[idim] = doTransformation1D(cdf[idim], *samples1d, bandwidths[idim], kern);

      // Update the kernel for the next dimension
      for (size_t isamples = 0; isamples < nsamples; isamples++) {
        xi = (cdf[idim] - samples1d->get(isamples)) / bandwidths[idim];
        kern[isamples] *=
            kde->getKernel().eval(xi);  // std::exp(-(xi * xi) / 2.); (bw*sqrt(2*PI)) cancels;
      }
    }

    // write them to the output
    pointsUniform.setRow(idata, unif);
  }

  return;
}

double OperationRosenblattTransformationKDE::doTransformation1D(double x, DataVector& samples1d,
                                                                double sigma, DataVector& kern) {
  // helper variables
  double cdfNormal = 0.0;
  double cdfConditionalized = 0.0;
  double xi = 0.;

  // compute dependent CDF over all kernels and build denominator
  double denom = 0.0;

  for (size_t isample = 0; isample < nsamples; isample++) {
    xi = (x - samples1d[isample]) / sigma;
    cdfNormal = kde->getKernel().cdf(xi);             // 0.5 + 0.5 * std::erf(xi / M_SQRT2);
    cdfConditionalized += kern[isample] * cdfNormal;  // (xx > xi(id,is));
    denom += kern[isample];

#ifdef DEBUG_ROSENBLATT
    std::cout << "i = " << isample << ": x_i = " << samples1d[isample] << " -> y_i = " << xi
              << std::endl;
    std::cout << "  cdfNormal = " << cdfNormal << ", denom = " << denom << std::endl;
#endif
  }

  // conditionalize the result
  cdfConditionalized /= denom;

  return cdfConditionalized;
}
}  // namespace datadriven
}  // namespace sgpp
