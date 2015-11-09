// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <map>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <random>
#include <algorithm>

#include "OperationRosenblattTransformationKDE.hpp"

#include <sgpp/globaldef.hpp>

//#ifdef _OPENMP
//#include <omp.h>
//#endif

//#define DEBUG_ROSENBLATT

using namespace SGPP::base;
using namespace SGPP::datadriven;

namespace SGPP {
  namespace datadriven {

    OperationRosenblattTransformationKDE::OperationRosenblattTransformationKDE(
      GaussianKDE& density) :
      kde(&density), bandwidths(density.getDim()), ndim(density.getDim()),
      nsamples(density.getNsamples()) {
      // get the bandwidth from the density for optimized
      density.getBandwidths(bandwidths);
    }

    OperationRosenblattTransformationKDE::~OperationRosenblattTransformationKDE() {
    }

    void OperationRosenblattTransformationKDE::doTransformation(
      DataMatrix& pointsCdf, DataMatrix& pointsUniform) {
      // Work arrays
      DataVector unif(ndim);
      DataVector cdf(ndim);
      DataVector kern(nsamples);
      DataVector* samples1d = nullptr;

      float_t xi = 0;

      for (size_t idata = 0; idata < pointsCdf.getNrows(); idata++) {
        kern.setAll(1.0);
        pointsCdf.getRow(idata, cdf);

        for (size_t idim = 0; idim < ndim; idim++) {
          // get samples in current dimension
          samples1d = kde->getSamples(idim);

          // transform the point in the current dimension
          unif[idim] = doTransformation1D(cdf[idim], *samples1d,
                                          bandwidths[idim], kern);

          // Update the kernel for the next dimension
          for (size_t isamples = 0; isamples < nsamples; isamples++) {
            xi = (cdf[idim] - samples1d->get(isamples)) / bandwidths[idim];
            kern[isamples] *= std::exp(-(xi * xi) / 2.); // (bw*sqrt(2*PI)) cancels;
          }
        }

        // write them to the output
        pointsUniform.setRow(idata, unif);
      }

      return;
    }

    void OperationRosenblattTransformationKDE::doShuffledTransformation(
      DataMatrix& pointsCdf, DataMatrix& pointsUniform) {
      // Work arrays
      DataVector unif(ndim);
      DataVector cdf(ndim);
      DataVector kern(nsamples);
      DataVector* samples1d = nullptr;

      float_t xi = 0;

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
          unif[idim] = doTransformation1D(cdf[idim], *samples1d,
                                          bandwidths[idim], kern);

          // Update the kernel for the next dimension
          for (size_t isamples = 0; isamples < nsamples; isamples++) {
            xi = (cdf[idim] - samples1d->get(isamples)) / bandwidths[idim];
            kern[isamples] *= std::exp(-(xi * xi) / 2.); // (bw*sqrt(2*PI)) cancels;
          }
        }

        // write them to the output
        pointsUniform.setRow(idata, unif);
      }

      return;
    }


    float_t OperationRosenblattTransformationKDE::doTransformation1D(float_t x,
        DataVector& samples1d, float_t sigma, DataVector& kern) {
      // helper variables
      float_t cdfNormal = 0.0;
      float_t cdfConditionalized = 0.0;
      float_t xi = 0.;

      // compute dependent CDF over all kernels and build denominator
      float_t denom = 0.0;

      for (size_t isample = 0; isample < nsamples; isample++) {
        xi = (x - samples1d[isample]) / sigma;
        cdfNormal = 0.5 + 0.5 * std::erf(xi / M_SQRT2);
        cdfConditionalized += kern[isample] * cdfNormal; // (xx > xi(id,is));
        denom += kern[isample];

#ifdef DEBUG_ROSENBLATT
        std::cout << "i = " << isample << ": x_i = " << samples1d[isample]
                  << " -> y_i = " << xi << std::endl;
        std::cout << "  cdfNormal = " << cdfNormal << ", denom = " << denom
                  << std::endl;
#endif
      }

      // conditionalize the result
      cdfConditionalized /= denom;

      return cdfConditionalized;
    }

  }
}
