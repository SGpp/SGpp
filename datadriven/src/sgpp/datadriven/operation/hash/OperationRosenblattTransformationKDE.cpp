/******************************************************************************
 * Copyright (C) 2009 Technische sample1diversitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de
#include <map>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <random>

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

//void OperationRosenblattTransformationKDE::doShuffeledTransformation(
//        DataMatrix& pointsUniform, DataMatrix& pointsCdf) {
//    // Work arrays
//    DataVector unif(ndim);
//    DataVector cdf(ndim);
//    DataVector kern(nsamples);
//    DataVector samples1d(nsamples);
//    size_t dim = 0;
//
//    // generate sample generator
//    std::random_device rd;
//    std::mt19937 generator(rd());
//    std::uniform_int_distribution<size_t> uniform_dist(0, ndim - 1);
//
//    // run over all samples
//    for (size_t isample = 0; isample < nsamples; isample++) {
//        kern.setAll(1.0);
//        pointsUniform.getRow(isample, unif);
//        // start in some random dimension
//        dim = uniform_dist(generator);
//        for (size_t idim = 0; idim < ndim; idim++) {
//            // transform the point in the current dimension
//            std::cout << dim << " ";
//
//            // get samples in current dimension
//            samples1d = samplesVec[idim];
//
//            // transform the point in the current dimension
//            unif[idim] = doTransformation1D(cdf[idim], samples1d,
//                    bandwidths[idim], kern);
//
//            dim = (dim + 1) % ndim;
//        }
//        // write them to the output
//        pointsCdf.setRow(isample, unif);
//        std::cout << std::endl;
//    }
//    return;
//}

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
