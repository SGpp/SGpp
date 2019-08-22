// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#include <sgpp/datadriven/operation/hash/simple/OperationCovariance.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalize.hpp>

#include <sgpp/base/operation/hash/OperationFirstMoment.hpp>
#include <sgpp/base/operation/hash/OperationSecondMoment.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/globaldef.hpp>
#include <vector>

namespace sgpp {
namespace datadriven {

base::DataMatrix* OperationCovariance::loadBounds(size_t numDims, base::DataMatrix* bounds,
                                                  size_t idim, size_t jdim) {
  base::DataMatrix* boundsdd = nullptr;
  if (bounds != nullptr) {
    boundsdd = new base::DataMatrix(numDims, 2);
    boundsdd->resize(numDims, 2);
    boundsdd->set(0, 0, bounds->get(idim, 0));
    boundsdd->set(0, 1, bounds->get(idim, 1));
    if (numDims > 1) {
      boundsdd->set(1, 0, bounds->get(jdim, 0));
      boundsdd->set(1, 1, bounds->get(jdim, 1));
    }
  }
  return boundsdd;
}

double OperationCovariance::mean(base::Grid& grid, base::DataVector& alpha,
                                 base::DataMatrix* bounds) {
  // compute the first moment given the boundaries
  std::unique_ptr<base::OperationFirstMoment> opFirstMoment(
      op_factory::createOperationFirstMoment(grid));
  return opFirstMoment->doQuadrature(alpha, bounds);
}

double OperationCovariance::variance(base::Grid& grid, base::DataVector& alpha,
                                     base::DataMatrix* bounds) {
  // compute the variance given the boundaries
  double firstMoment = mean(grid, alpha, bounds);
  std::unique_ptr<base::OperationSecondMoment> opSecondMoment(
      op_factory::createOperationSecondMoment(grid));
  double secondMoment = opSecondMoment->doQuadrature(alpha, bounds);
  // use Steiners translation theorem to compute the variance
  return secondMoment - firstMoment * firstMoment;
}

void OperationCovariance::doQuadrature(base::DataVector& alpha, base::DataMatrix& cov,
                                       base::DataMatrix* bounds) {
  size_t ndim = grid.getDimension();

  if ((cov.getNrows() != ndim) || (cov.getNcols() != ndim)) {
    // covariance matrix has wrong size -> resize
    cov.resize(ndim, ndim);
  }

  // prepare covariance marix
  cov.setAll(0.0);

  // generate 1d densities and compute means and variances
  base::DataVector means(ndim);
  base::DataVector variances(ndim);

  std::unique_ptr<datadriven::OperationDensityMargTo1D> opMarg(
      op_factory::createOperationDensityMargTo1D(grid));

  base::Grid* marginalizedGrid = nullptr;
  base::DataVector* marginalizedAlpha = new base::DataVector(0);

  for (size_t idim = 0; idim < ndim; idim++) {
    // marginalize and normalize
    opMarg->margToDimX(&alpha, marginalizedGrid, marginalizedAlpha, idim);
    // load bounding box and compute moments
    base::DataMatrix* boundsdd = loadBounds(1, bounds, idim);
    means[idim] = mean(*marginalizedGrid, *marginalizedAlpha, boundsdd);
    variances[idim] = variance(*marginalizedGrid, *marginalizedAlpha, boundsdd);
    delete marginalizedGrid;
    if (boundsdd != nullptr) {
      delete boundsdd;
    }
  }

  // helper variables
  std::vector<size_t> mdims(2);
  double covij = 0.0;
  for (size_t idim = 0; idim < ndim; idim++) {
    // diagonal is equal to the variance of the marginalized densities
    cov.set(idim, idim, variances[idim]);
    for (size_t jdim = idim + 1; jdim < ndim; jdim++) {
      // marginalize the density
      mdims[0] = idim;
      mdims[1] = jdim;
      opMarg->margToDimXs(&alpha, marginalizedGrid, marginalizedAlpha, mdims);
      // -----------------------------------------------------
      // compute the covariance of Cov(X_i, X_j)
      base::DataMatrix* boundsdd = loadBounds(2, bounds, idim, jdim);
      covij = mean(*marginalizedGrid, *marginalizedAlpha, boundsdd) - means[idim] * means[jdim];
      cov.set(idim, jdim, covij);
      cov.set(jdim, idim, covij);
      // -----------------------------------------------------
      delete marginalizedGrid;
      if (boundsdd != nullptr) {
        delete boundsdd;
      }
    }
  }

  delete marginalizedAlpha;
}

}  // namespace datadriven
}  // namespace sgpp
