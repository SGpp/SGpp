/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de
// some defines for the following algorithm
#include <map>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <random>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include "OperationInverseRosenblattTransformationKDE.hpp"
#include "OperationRosenblattTransformationKDE.hpp"

#include <sgpp/globaldef.hpp>

//#ifdef _OPENMP
//#include <omp.h>
//#endif

//#define DEBUG_INVERSE_ROSENBLATT

using namespace SGPP::base;
using namespace SGPP::datadriven;

namespace SGPP {
  namespace datadriven {

    // ----------------------------------------------------------------------------

    OperationInverseRosenblattTransformationKDE::OperationInverseRosenblattTransformationKDE(
      GaussianKDE& kde, float_t sigmaFactor, float_t inversionEpsilon) :
      kde(&kde), bandwidths(kde.getDim()), xlimits(2, kde.getDim()), ylimits(
        2, kde.getDim()), inversionEpsilon(inversionEpsilon) {
      // estimate bandwidth and copy result to local sigma
      kde.getBandwidths(bandwidths);

      ndim = kde.getDim();
      nsamples = kde.getNsamples();

      // define search limits based on the data and the estimated bandwidths
      // numerical instabilities expected for sigmaFactor > 9
      recalcLimits(sigmaFactor);
    }

    OperationInverseRosenblattTransformationKDE::~OperationInverseRosenblattTransformationKDE() {
    }

    // ----------------------------------------------------------------------------

    float_t OperationInverseRosenblattTransformationKDE::getMaxInversionError() {
      return inversionEpsilon;
    }

    void OperationInverseRosenblattTransformationKDE::recalcLimits(
      float_t sigmaFactor) {
      DataVector* samples1d = nullptr;

      xlimits.resize(2, ndim);
      ylimits.resize(2, ndim);

      DataVector xlimits_1d(2);
      DataVector ylimits_1d(2);

      DataVector kern(nsamples);
      kern.setAll(1.0f);

      OperationRosenblattTransformationKDE* opRosen =
        op_factory::createOperationRosenblattTransformationKDE(*kde);

      // get minimum and maximum of the data in every dimension
      for (size_t idim = 0; idim < ndim; idim++) {
        samples1d = kde->getSamples(idim);

        // evaluate the confidence interval for given alpha value
        xlimits_1d[0] = samples1d->min() - sigmaFactor * bandwidths[idim];
        xlimits_1d[1] = samples1d->max() + sigmaFactor * bandwidths[idim];

        // transform the limits to obtain the limits in [0, 1] for inversion
        ylimits_1d[0] = opRosen->doTransformation1D(xlimits_1d[0], *samples1d,
                        bandwidths[idim], kern);
        ylimits_1d[1] = opRosen->doTransformation1D(xlimits_1d[1], *samples1d,
                        bandwidths[idim], kern);

        // store the obtained limits
        xlimits.setColumn(idim, xlimits_1d);
        ylimits.setColumn(idim, ylimits_1d);

#ifdef DEBUG_INVERSE_ROSENBLATT
        std::cout << "d=" << idim << ":" << std::endl;
        std::cout << " data in [" << samples1d->min() << ", "
                  << samples1d->max() << "], length = " << samples1d->getSize()
                  << std::endl;
        std::cout << " sigma = " << bandwidths[idim] << std::endl;
        std::cout << " x    in [" << xlimits_1d[0] << ", " << xlimits_1d[1]
                  << "]" << std::endl;
        std::cout << " y    in [" << ylimits_1d[0] << ", " << ylimits_1d[1]
                  << "]" << std::endl;
#endif
      }

      delete opRosen;

      return;
    }

    void OperationInverseRosenblattTransformationKDE::doTransformation(
      DataMatrix& pointsUniform, DataMatrix& pointsCdf) {
      // Work arrays
      DataVector unif(ndim);
      DataVector cdf(ndim);
      DataVector kern(nsamples);
      DataVector* samples1d = nullptr;

      float_t xi = 0;

      for (size_t idata = 0; idata < pointsUniform.getNrows(); idata++) {
        kern.setAll(1.0);
        pointsUniform.getRow(idata, unif);

        for (size_t idim = 0; idim < ndim; idim++) {
          // get samples in current dimension
          samples1d = kde->getSamples(idim);

          // transform the point in the current dimension
          cdf[idim] = doTransformation1D(unif[idim], *samples1d,
                                         bandwidths[idim], xlimits.get(0, idim),
                                         xlimits.get(1, idim), ylimits.get(0, idim),
                                         ylimits.get(1, idim), kern);

          // Update the kernel for the next dimension
          for (size_t isamples = 0; isamples < nsamples; isamples++) {
            xi = (cdf[idim] - samples1d->get(isamples)) / bandwidths[idim];
            kern[isamples] *= exp(-(xi * xi) / 2.); // (bw*sqrt(2*PI)) cancels;
          }
        }

        // write them to the output
        pointsCdf.setRow(idata, cdf);
      }

      return;
    }

    //void OperationInverseRosenblattTransformationKDE::doShuffeledTransformation(
    //        base::DataMatrix* pointsUniform, base::DataMatrix* pointsCdf) {
    //    size_t ndim = this->trainingData->getNcols();
    //    size_t nsamples = this->trainingData->getNrows();
    //
    //    // Work arrays
    //    DataVector unif(ndim);
    //    DataVector kern(nsamples);
    //    float_t uniform1d = 0.0;
    //    float_t cdf1d = 0.0;
    //    size_t dim = 0;
    //
    //    // generate sample generator
    //    std::random_device rd;
    //    std::mt19937 generator(rd());
    //    std::uniform_int_distribution<int> uniform_dist(0, ndim - 1);
    //
    //    // run over all samples
    //    for (size_t i = 0; i < nsamples; i++) {
    //        kern.setAll(1.0);
    //        pointsUniform->getRow(i, unif);
    //        // start in some random dimension
    //        dim = uniform_dist(generator);
    //        for (size_t j = 0; j < ndim; j++) {
    //            // transform the point in the current dimension
    //            std::cout << dim << " ";
    //            unif[dim] = doTransformation1D(unif[dim], kern, dim);
    //            dim = (dim + 1) % ndim;
    //        }
    //        // write them to the outpu
    //        pointsCdf->setRow(i, unif);
    //        std::cout << std::endl;
    //    }
    //    return;
    //}

    //void OperationInverseRosenblattTransformationKDE::doTransformation_start_dimX(
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
    //void OperationInverseRosenblattTransformationKDE::doTransformation_in_next_dim(
    //        Grid* g_in, DataVector* a_in, size_t op_dim,
    //        DataVector* cdfs1d, DataVector* coords1d,
    //        size_t& curr_dim) {
    //
    //    size_t dims = cdfs1d->getSize();  // total dimensions
    //
    //    /* Step 1: do conditional in current dim */
    //    Grid* g_out = NULL;
    //    DataVector* a_out = new DataVector(1);
    //    OperationDensityConditional* cond =
    //            op_factory::createOperationDensityConditional(*g_in);
    //    cond->doConditional(*a_in, g_out, *a_out, static_cast<unsigned int>(op_dim),
    //            coords1d->get(curr_dim));
    //    delete cond;
    //
    //    // move on to next dim
    //    curr_dim = (curr_dim + 1) % dims;
    //    op_dim = (op_dim + 1) % g_out->getStorage()->dim();
    //
    //    /* Step 2: draw a sample in next dim */
    //    float_t x = 0;
    //    if (g_out->getStorage()->dim() > 1) {
    //
    //        // Marginalize to next dimension
    //        Grid* g1d = NULL;
    //        DataVector* a1d = NULL;
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
    //    if (g_out->getStorage()->dim() > 1)
    //        doTransformation_in_next_dim(g_out, a_out, op_dim, cdfs1d, coords1d,
    //                curr_dim);
    //
    //    delete g_out;
    //    delete a_out;
    //
    //    return;
    //}

    float_t OperationInverseRosenblattTransformationKDE::doTransformation1D(
      float_t coord1d, DataVector& samples1d, float_t sigma, float_t xlower,
      float_t xupper, float_t ylower, float_t yupper, DataVector& kern) {
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
      float_t denom = 0.0;

      for (size_t isamples = 0; isamples < nsamples; isamples++) {
        denom += kern[isamples];
      }

      // stores the result
      float_t xres = 0.;

      // apply bisection method to obtain a starting posize_t for newton's method
      float_t xerr = 0.0;
      float_t xacc = 1e-1;

      float_t xBisection = 0.0;
      xerr = bisection(coord1d, xBisection, xlower, xupper, samples1d, sigma,
                       kern, denom, xacc);
      // use this as starting point for next solver
      float_t xNext = xBisection;

      // use this as starting point for next solver (newton) with higher accuracy
      xerr = newton(coord1d, xNext, samples1d, sigma, kern, denom,
                    inversionEpsilon);
      //    xerr = halley(coord1d, xNext, dim, kern, denom, xacc);

      if (std::isnan(xNext) || xerr > inversionEpsilon) {
        // newton is not converged -> run bisection again with higher accuracy
        xerr = bisection(coord1d, xBisection, xlower, xupper, samples1d, sigma,
                         kern, denom, inversionEpsilon);
        xres = xBisection;

        if (xerr > xacc) {
          throw algorithm_exception(
            "Error: inversion with Rosenblatt is not converged. Search interval for root is possibly too small.");
        }
      } else {
        xres = xNext;
      }

      return xres;
    } // end of compute_1D_cdf()

    float_t OperationInverseRosenblattTransformationKDE::bisection(float_t y,
        float_t& x, float_t& xlower, float_t& xupper,
        base::DataVector& samples1d, float_t sigma, DataVector& kern,
        float_t denom, float_t xacc, size_t maxIterations) {
      // iteration counter
      size_t ii = 0;
      size_t ns = samples1d.getSize();

      // use Bisection as starting posize_t for newton's method
      // Bisection looks in the size_terval [x1,x2]
      float_t cdfNormal = 0.0;
      float_t cdfConditionalized = 0.0;
      float_t xi = 0.0;
      float_t xerr = 0.;

      // Bisection loop = search zeros in order to invert CDF = erf
      do {
        cdfConditionalized = 0;

        // compute dependent CDF over all kernels
        for (size_t is = 0; is < ns; is++) {
          xi = (x - samples1d[is]) / sigma;
          cdfNormal = 0.5 + 0.5 * erf(xi / M_SQRT2);
          cdfConditionalized += kern[is] * cdfNormal; // (xx > xi(id,is));
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
        std::cout << "Warning: bisection not converged for x " << " in ["
                  << xlower << ", " << xupper << "], err=" << xerr << " > "
                  << xacc << ", iterations=" << ii << std::endl;
      } else {
        std::cout << "bisection: x = " << x << ", err=" << xerr << " < " << xacc
                  << ", iterations=" << ii << "/" << maxIterations << std::endl;
      }

#endif

      return xerr;
    }

    float_t OperationInverseRosenblattTransformationKDE::newton(float_t y,
        float_t& x, base::DataVector& samples1d, float_t sigma,
        DataVector& kern, float_t denom, float_t xacc, size_t maxIterations) {
      // newton method
      size_t ii = 0;
      size_t ns = samples1d.getSize();

      // helper variables
      float_t cdfNormal = 0.0;
      float_t pdfNormal = 0.;

      // conditionalized cdf function value
      float_t fx = 0.0;
      // conditionalized pdf function value = derivative of cdf
      float_t fdx = 0.0;
      // normalized normal
      float_t xi = 0.0;
      // stores old x value and current error
      float_t xold = 0;
      float_t xerr = 0.0;

      do {
        // store old x-value
        xold = x;

        fx = 0;
        fdx = 0;

        // compute dependent CDF over all kernels
        for (size_t is = 0; is < ns; is++) {
          xi = (x - samples1d[is]) / sigma;
          pdfNormal = std::exp(-(xi * xi) / 2.);
          cdfNormal = 0.5 + 0.5 * erf(xi / M_SQRT2);

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
        std::cout << "Warning: newton not converged for x " << ", err=" << xerr
                  << " > " << xacc << ", iterations=" << ii << std::endl;
      } else {
        std::cout << "newton   : x = " << x << ", err=" << xerr << " < " << xacc
                  << ", iterations=" << ii << "/" << maxIterations << std::endl;
      }

#endif
      return xerr;
    }

    float_t OperationInverseRosenblattTransformationKDE::halley(float_t y,
        float_t& x, base::DataVector& samples1d, float_t sigma,
        DataVector& kern, float_t denom, float_t xacc, size_t maxIterations) {
      size_t ii = 0;
      size_t ns = samples1d.getSize();

      float_t cdfNormal = 0.0;
      float_t pdfNormal = 0.;
      float_t derivativePdfNormal = 0;

      // conditionalized cdf function value
      float_t fx = 0.0;
      // conditionalized pdf function value = derivative of cdf
      float_t fdx = 0.0;
      // conditionalized derivative of pdf function value
      float_t fddx = 0.0;

      // normalized normal
      float_t xi = 0.0;

      // stores old x value and current error
      float_t xold = 0;
      float_t xerr = 0.0;

      do {
        // store old x-value
        xold = x;

        fx = 0;
        fdx = 0;
        fddx = 0.0;

        // compute dependent CDF over all kernels
        for (size_t is = 0; is < ns; is++) {
          xi = (x - samples1d[is]) / sigma;
          cdfNormal = 0.5 + 0.5 * erf(xi / M_SQRT2);
          pdfNormal = exp(-(xi * xi) / 2.);
          derivativePdfNormal = xi / sigma * pdfNormal;

          fx += kern[is] * cdfNormal;
          fdx += kern[is] * pdfNormal;
          fddx += kern[is] * derivativePdfNormal;
        }

        fdx *= M_1_SQRT2PI / sigma;
        fddx *= M_1_SQRT2PI / sigma;

        fx = y - fx / denom;
        fdx = -fdx / denom;
        fddx = -fddx / denom;

        // apply halleys method
        x -= (2 * fx * fdx) / (2 * fdx * fdx - fx * fddx);
        xerr = fabs(x - xold);
        ii++;
      } while ((xerr > xacc) && (ii < maxIterations));

      // end of halleys method

#ifdef DEBUG_INVERSE_ROSENBLATT

      if (ii >= maxIterations || std::isnan(x) || xerr > xacc) {
        std::cout << "Warning: halley not converged, err=" << xerr << " > "
                  << xacc << ", iterations=" << ii << std::endl;
      } else {
        std::cout << "halley: err=" << xerr << " < " << xacc << ", iterations="
                  << ii << "/" << maxIterations << std::endl;
      }

#endif
      return xerr;
    }

  }
}
