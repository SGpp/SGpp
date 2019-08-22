/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de
// some defines for the following algorithm
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationKDE.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationKDE.hpp>

#include <sgpp/globaldef.hpp>

#include <map>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// #define DEBUG_INVERSE_ROSENBLATT

namespace sgpp {
namespace datadriven {

// ----------------------------------------------------------------------------

OperationInverseRosenblattTransformationKDE::OperationInverseRosenblattTransformationKDE(
    KernelDensityEstimator& kde, double sigmaFactor, double inversionEpsilon, std::uint64_t seed)
    : kde(&kde),
      bandwidths(kde.getDim()),
      xlimits(2, kde.getDim()),
      ylimits(2, kde.getDim()),
      inversionEpsilon(inversionEpsilon),
      rng(seed) {
  // estimate bandwidth and copy result to local sigma
  kde.getBandwidths(bandwidths);

  ndim = kde.getDim();
  nsamples = kde.getNsamples();

  // define search limits based on the data and the estimated bandwidths
  // numerical instabilities expected for sigmaFactor > 9
  recalcLimits(sigmaFactor);
}

OperationInverseRosenblattTransformationKDE::~OperationInverseRosenblattTransformationKDE() {}

// ----------------------------------------------------------------------------

double OperationInverseRosenblattTransformationKDE::getMaxInversionError() {
  return inversionEpsilon;
}

void OperationInverseRosenblattTransformationKDE::recalcLimits(double sigmaFactor) {
  std::shared_ptr<base::DataVector> samples1d;

  xlimits.resize(2, ndim);
  ylimits.resize(2, ndim);

  base::DataVector xlimits_1d(2);
  base::DataVector ylimits_1d(2);

  base::DataVector kern(nsamples);
  kern.setAll(1.0);

  std::unique_ptr<OperationRosenblattTransformationKDE> opRosen(
      op_factory::createOperationRosenblattTransformationKDE(*kde));

  // get minimum and maximum of the data in every dimension
  for (size_t idim = 0; idim < ndim; idim++) {
    samples1d = kde->getSamples(idim);

    // evaluate the confidence interval for given alpha value
    xlimits_1d[0] = samples1d->min() - sigmaFactor * bandwidths[idim];
    xlimits_1d[1] = samples1d->max() + sigmaFactor * bandwidths[idim];

    // transform the limits to obtain the limits in [0, 1] for inversion
    ylimits_1d[0] = opRosen->doTransformation1D(xlimits_1d[0], *samples1d, bandwidths[idim], kern);
    ylimits_1d[1] = opRosen->doTransformation1D(xlimits_1d[1], *samples1d, bandwidths[idim], kern);

    // store the obtained limits
    xlimits.setColumn(idim, xlimits_1d);
    ylimits.setColumn(idim, ylimits_1d);

#ifdef DEBUG_INVERSE_ROSENBLATT
    std::cout << "d=" << idim << ":" << std::endl;
    std::cout << " data in [" << samples1d->min() << ", " << samples1d->max()
              << "], length = " << samples1d->getSize() << std::endl;
    std::cout << " sigma = " << bandwidths[idim] << std::endl;
    std::cout << " x    in [" << xlimits_1d[0] << ", " << xlimits_1d[1] << "]" << std::endl;
    std::cout << " y    in [" << ylimits_1d[0] << ", " << ylimits_1d[1] << "]" << std::endl;
#endif
  }

  return;
}

void OperationInverseRosenblattTransformationKDE::doTransformation(base::DataMatrix& pointsUniform,
                                                                   base::DataMatrix& pointsCdf) {
#pragma omp parallel
  {
#pragma omp for schedule(dynamic)
    for (size_t idata = 0; idata < pointsUniform.getNrows(); idata++) {
      // Work arrays
      base::DataVector unif(ndim);
      base::DataVector cdf(ndim);
      base::DataVector kern(nsamples);
      std::shared_ptr<base::DataVector> samples1d;

      kern.setAll(1.0);
      pointsUniform.getRow(idata, unif);

      for (size_t idim = 0; idim < ndim; idim++) {
        // get samples in current dimension
        samples1d = kde->getSamples(idim);

        // transform the point in the current dimension
        cdf[idim] = doTransformation1D(unif[idim], *samples1d, bandwidths[idim],
                                       xlimits.get(0, idim), xlimits.get(1, idim),
                                       ylimits.get(0, idim), ylimits.get(1, idim), kern);

        // Update the kernel for the next dimension
        for (size_t isamples = 0; isamples < nsamples; isamples++) {
          double xi = (cdf[idim] - samples1d->get(isamples)) / bandwidths[idim];
          kern[isamples] *=
              kde->getKernel().eval(xi);  // exp(-(xi * xi) / 2.);  (bw*sqrt(2*PI)) cancels;
        }
      }

      // write them to the output
      pointsCdf.setRow(idata, cdf);
    }
  }
  return;
}

void OperationInverseRosenblattTransformationKDE::doShuffledTransformation(
    base::DataMatrix& pointsUniform, base::DataMatrix& pointsCdf) {
  // 1. compute permuations for each sample
  size_t num_dims = this->kde->getDim();
  size_t num_samples = pointsCdf.getNrows();
  std::vector<std::vector<size_t>> permutations(num_samples);

  std::vector<size_t> startindices(num_samples);
  // change the starting dimension when the bucket_size is arrived
  // this distributes the error in the projection uniformly to all
  // dimensions and make it therefore stable
  size_t dim_start = 0;
  size_t bucket_size = num_samples / num_dims + 1;
  for (size_t i = 0; i < num_samples; i++) {
    if (((i + 1) % bucket_size) == 0 && (i + 1) < num_samples) {
      ++dim_start;
    }
    permutations[i].resize(ndim);
    for (size_t idim = 0; idim < ndim; idim++) {
      permutations[i][idim] = (dim_start + idim) % num_dims;
    }
  }

// apply the inverse rosenblatt transformation
#pragma omp parallel
  {
#pragma omp for schedule(dynamic)
    for (size_t idata = 0; idata < pointsCdf.getNrows(); idata++) {
      // Work arrays
      base::DataVector unif(ndim);
      base::DataVector cdf(ndim);
      base::DataVector kern(nsamples);
      std::shared_ptr<base::DataVector> samples1d;

      kern.setAll(1.0);
      pointsUniform.getRow(idata, unif);

      for (size_t i = 0; i < ndim; i++) {
        size_t idim = permutations[idata][i];
        // get samples in current dimension
        samples1d = kde->getSamples(idim);

        // transform the point in the current dimension
        cdf[idim] = doTransformation1D(unif[idim], *samples1d, bandwidths[idim],
                                       xlimits.get(0, idim), xlimits.get(1, idim),
                                       ylimits.get(0, idim), ylimits.get(1, idim), kern);

        // Update the kernel for the next dimension
        for (size_t isamples = 0; isamples < nsamples; isamples++) {
          double xi = (cdf[idim] - samples1d->get(isamples)) / bandwidths[idim];
          kern[isamples] *=
              kde->getKernel().eval(xi);  // exp(-(xi * xi) / 2.); (bw*sqrt(2*PI)) cancels;
        }
      }

      // write them to the output
      pointsCdf.setRow(idata, cdf);
    }
  }
  return;
}

// void OperationInverseRosenblattTransformationKDE::doTransformation_start_dimX(
//        Grid* g_in, DataVector* a_in, size_t dim_start,
//        DataVector* cdfs1d, DataVector* coords1d) {
//
//    size_t dims = coords1d->getSize(); // total dimensions
//
//    if ((dims > 1) && (dim_start <= dims - 1)) {
//        size_t curr_dim = dim_start;
//        doTransformation_in_next_dim(g_in, a_in, dim_start, cdfs1d, coords1d,
//                curr_dim);
//    } else if (dims == 1) {
//        throw operation_exception(
//                "Error: # of dimensions = 1. No operation needed!");
//    } else {
//        throw operation_exception(
//                "Error: dimension out of range. Operation aborted!");
//    }
//
//    return;
//}
//
// void OperationInverseRosenblattTransformationKDE::doTransformation_in_next_dim(
//        Grid* g_in, DataVector* a_in, size_t op_dim,
//        DataVector* cdfs1d, DataVector* coords1d,
//        size_t& curr_dim) {
//
//    size_t dims = cdfs1d->getSize();  // total dimensions
//
//    /* Step 1: do conditional in current dim */
//    Grid* g_out = nullptr;
//    DataVector* a_out = new DataVector(1);
//    OperationDensityConditional* cond =
//            op_factory::createOperationDensityConditional(*g_in);
//    cond->doConditional(*a_in, g_out, *a_out, static_cast<unsigned int>(op_dim),
//            coords1d->get(curr_dim));
//    delete cond;
//
//    // move on to next dim
//    curr_dim = (curr_dim + 1) % dims;
//    op_dim = (op_dim + 1) % g_out->getDimension();
//
//    /* Step 2: draw a sample in next dim */
//    double x = 0;
//    if (g_out->getDimension() > 1) {
//
//        // Marginalize to next dimension
//        Grid* g1d = nullptr;
//        DataVector* a1d = nullptr;
//        OperationDensityMargTo1D* marg1d =
//                op_factory::createOperationDensityMargTo1D(*g_out);
//        marg1d->margToDimX(a_out, g1d, a1d, op_dim);
//        delete marg1d;
//
//        // Draw a sample in next dimension
//        x = doTransformation1D(g1d, a1d, cdfs1d->get(curr_dim));
//        delete g1d;
//        delete a1d;
//
//    } else {
//        // skip Marginalize, directly draw a sample in next dimension
//        x = doTransformation1D(g_out, a_out, cdfs1d->get(curr_dim));
//    }
//
//    /* Step 3: copy sample to output */
//    coords1d->set(curr_dim, x);
//
//    /* Step 4: sample in next dimension */
//    if (g_out->getDimension() > 1)
//        doTransformation_in_next_dim(g_out, a_out, op_dim, cdfs1d, coords1d,
//                curr_dim);
//
//    delete g_out;
//    delete a_out;
//
//    return;
//}

double OperationInverseRosenblattTransformationKDE::doTransformation1D(
    double coord1d, base::DataVector& samples1d, double sigma, double xlower, double xupper,
    double ylower, double yupper, base::DataVector& kern) {
  // Cure against extremes
  if (coord1d <= ylower) {
    return xlower;
  }

  if (coord1d >= yupper) {
    return xupper;
  }

  // get the number of samples
  size_t nsamples = kern.getSize();

  // Build denominator
  double denom = 0.0;

  for (size_t isamples = 0; isamples < nsamples; isamples++) {
    denom += kern[isamples];
  }

  // stores the result
  double xres = 0.;

  // apply bisection method to obtain a starting posize_t for newton's method
  double xerr = 0.0;
  double xacc = 1e-1;

  double xBisection = 0.0;
  xerr = bisection(coord1d, xBisection, xlower, xupper, samples1d, sigma, kern, denom, xacc);
  // use this as starting point for next solver
  double xNext = xBisection;

  // use this as starting point for next solver (newton) with higher accuracy
  xerr = newton(coord1d, xNext, samples1d, sigma, kern, denom, inversionEpsilon);

  if (std::isnan(xNext) || xerr > inversionEpsilon) {
    // newton is not converged -> run bisection again with higher accuracy
    xerr = bisection(coord1d, xBisection, xlower, xupper, samples1d, sigma, kern, denom,
                     inversionEpsilon);
    xres = xBisection;

    if (xerr > xacc) {
      throw base::algorithm_exception(
          "Error: inversion with Rosenblatt is not converged. Search interval for root is possibly "
          "too small.");
    }
  } else {
    xres = xNext;
  }

  return xres;
}  // end of compute_1D_cdf()

double OperationInverseRosenblattTransformationKDE::bisection(
    double y, double& x, double& xlower, double& xupper, base::DataVector& samples1d, double sigma,
    base::DataVector& kern, double denom, double xacc, size_t maxIterations) {
  // iteration counter
  size_t ii = 0;
  size_t ns = samples1d.getSize();

  // use Bisection as starting posize_t for newton's method
  // Bisection looks in the size_terval [x1,x2]
  double cdfNormal = 0.0;
  double cdfConditionalized = 0.0;
  double xi = 0.0;
  double xerr = 0.;

  // Bisection loop = search zeros in order to invert CDF = erf
  do {
    cdfConditionalized = 0;

    // compute dependent CDF over all kernels
    for (size_t is = 0; is < ns; is++) {
      xi = (x - samples1d[is]) / sigma;
      cdfNormal = kde->getKernel().cdf(xi);        // 0.5 + 0.5 * erf(xi / M_SQRT2);
      cdfConditionalized += kern[is] * cdfNormal;  // (xx > xi(id,is));
    }

    // divide result by denominator
    cdfConditionalized /= denom;

    // adjust domain for possible zero
    if (y < cdfConditionalized) {
      xupper = x;
    } else {
      xlower = x;
    }

    // select new center posize_t
    x = (xlower + xupper) / 2.0;
    xerr = xupper - xlower;
    ii++;
  } while ((xerr > xacc) && (ii < maxIterations));

#ifdef DEBUG_INVERSE_ROSENBLATT

  if (ii >= maxIterations || std::isnan(x) || xerr > xacc) {
    std::cout << "Warning: bisection not converged for x "
              << " in [" << xlower << ", " << xupper << "], err=" << xerr << " > " << xacc
              << ", iterations=" << ii << std::endl;
  } else {
    std::cout << "bisection: x = " << x << ", err=" << xerr << " < " << xacc
              << ", iterations=" << ii << "/" << maxIterations << std::endl;
  }

#endif

  return xerr;
}

double OperationInverseRosenblattTransformationKDE::newton(double y, double& x,
                                                           base::DataVector& samples1d,
                                                           double sigma, base::DataVector& kern,
                                                           double denom, double xacc,
                                                           size_t maxIterations) {
  // newton method
  size_t ii = 0;
  size_t ns = samples1d.getSize();

  // helper variables
  double cdfNormal = 0.0;
  double pdfNormal = 0.;

  // conditionalized cdf function value
  double fx = 0.0;
  // conditionalized pdf function value = derivative of cdf
  double fdx = 0.0;
  // normalized normal
  double xi = 0.0;
  // stores old x value and current error
  double xold = 0;
  double xerr = 0.0;

  do {
    // store old x-value
    xold = x;

    fx = 0;
    fdx = 0;

    // compute dependent CDF over all kernels
    for (size_t is = 0; is < ns; is++) {
      xi = (x - samples1d[is]) / sigma;
      pdfNormal = kde->getKernel().eval(xi);  // std::exp(-(xi * xi) / 2.);
      cdfNormal = kde->getKernel().cdf(xi);   // 0.5 + 0.5 * erf(xi / M_SQRT2);

      fx += kern[is] * cdfNormal;
      fdx += kern[is] * pdfNormal;
    }

    fdx *= M_1_SQRT2PI / sigma;

    // apply newtons method
    x += y * (denom / fdx) - fx / fdx;
    xerr = fabs(x - xold);
    ii++;
  } while ((xerr > xacc) && (ii < maxIterations));

// end of newton method

#ifdef DEBUG_INVERSE_ROSENBLATT

  if (ii >= maxIterations || std::isnan(x) || xerr > xacc) {
    std::cout << "Warning: newton not converged for x "
              << ", err=" << xerr << " > " << xacc << ", iterations=" << ii << std::endl;
  } else {
    std::cout << "newton   : x = " << x << ", err=" << xerr << " < " << xacc
              << ", iterations=" << ii << "/" << maxIterations << std::endl;
  }

#endif
  return xerr;
}
}  // namespace datadriven
}  // namespace sgpp
