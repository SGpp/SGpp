// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationModPolyClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityConditional.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensitySampling1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DModPolyClenshawCurtis.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <map>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sgpp {
namespace datadriven {

void OperationRosenblattTransformationModPolyClenshawCurtis::
     doTransformation(base::DataVector* alpha, base::DataMatrix* points,
                      base::DataMatrix* pointscdf) {
  size_t num_dims = this->grid->getDimension();

  // 1. marginalize to all possible start dimensions
  std::vector<base::Grid*> grids1d(num_dims);
  std::vector<base::DataVector*> alphas1d(num_dims);
  std::unique_ptr<OperationDensityMargTo1D> marg1d(
      op_factory::createOperationDensityMargTo1D(*this->grid));
  for (size_t idim = 0; idim < num_dims; idim++) {
    marg1d->margToDimX(alpha, grids1d[idim], alphas1d[idim], idim);
  }

  // 2. compute the start dimension for each sample
  size_t num_samples = pointscdf->getNrows();
  std::vector<size_t> startindices(num_samples);
  // change the starting dimension when the bucket_size is arrived
  // this distributes the error in the projection uniformly to all
  // dimensions and make it therefore stable
  size_t dim_start = 0;
  size_t bucket_size = num_samples / num_dims + 1;
  for (size_t i = 0; i < num_samples; i++) {
    if (((i + 1) % bucket_size) == 0 && (i + 1) < pointscdf->getNrows()) {
      ++dim_start;
    }
    startindices[i] = dim_start;
  }

// 3. for every sample do...
#pragma omp parallel
  {
#pragma omp for schedule(dynamic)
    for (size_t i = 0; i < points->getNrows(); i++) {
      // transform the point in the current dimension
      size_t idim = startindices[i];
      double y = doTransformation1D(grids1d[idim], alphas1d[idim], points->get(i, idim));
      // and write it to the output
      pointscdf->set(i, idim, y);

      // prepare the next dimensions -> read samples
      base::DataVector cdfs1d(num_dims);
      base::DataVector coords1d(num_dims);

      points->getRow(i, coords1d);
      pointscdf->getRow(i, cdfs1d);
      doTransformation_start_dimX(this->grid, alpha, idim, &coords1d, &cdfs1d);
      pointscdf->setRow(i, cdfs1d);
    }
  }

  // cleanup
  for (size_t idim = 0; idim < num_dims; idim++) {
    delete grids1d[idim];
    delete alphas1d[idim];
  }
}

void OperationRosenblattTransformationModPolyClenshawCurtis::
     doTransformation(base::DataVector* alpha, base::DataMatrix* points,
                      base::DataMatrix* pointscdf, size_t dim_start) {
  // 1. marginalize to dim_start
  base::Grid* g1d = nullptr;
  base::DataVector* a1d = nullptr;
  std::unique_ptr<OperationDensityMargTo1D> marg1d(
      op_factory::createOperationDensityMargTo1D(*this->grid));
  marg1d->margToDimX(alpha, g1d, a1d, dim_start);

#pragma omp parallel
  {
#pragma omp for schedule(dynamic)

    for (size_t i = 0; i < points->getNrows(); i++) {
      base::DataVector coords1d(points->getNcols());
      base::DataVector cdfs1d(points->getNcols());
      // 2. 1D transformation on dim_start
      double y = doTransformation1D(g1d, a1d, points->get(i, dim_start));
      pointscdf->set(i, dim_start, y);
      // 3. for every missing dimension do...
      points->getRow(i, coords1d);
      pointscdf->getRow(i, cdfs1d);
      doTransformation_start_dimX(this->grid, alpha, dim_start, &coords1d, &cdfs1d);
      pointscdf->setRow(i, cdfs1d);
    }
  }

  delete g1d;
  delete a1d;
}

void OperationRosenblattTransformationModPolyClenshawCurtis::doTransformation_start_dimX(
    base::Grid* g_in, base::DataVector* a_in, size_t dim_start, base::DataVector* coords1d,
    base::DataVector* cdfs1d) {
  size_t dims = coords1d->getSize();  // total dimensions

  if ((dims > 1) && (dim_start <= dims - 1)) {
    size_t curr_dim = dim_start;
    doTransformation_in_next_dim(g_in, a_in, dim_start, coords1d, cdfs1d, curr_dim);
  } else if (dims == 1) {
    throw base::operation_exception("Error: # of dimensions = 1. No operation needed!");
  } else {
    throw base::operation_exception("Error: dimension out of range. Operation aborted!");
  }

  return;
}

void OperationRosenblattTransformationModPolyClenshawCurtis::doTransformation_in_next_dim(
    base::Grid* g_in, base::DataVector* a_in, size_t op_dim, base::DataVector* coords1d,
    base::DataVector* cdfs1d, size_t& curr_dim) {
  size_t dims = coords1d->getSize();  // total dimensions
  // unsigned int op_dim = curr_dim + 1;  // (curr_dim < dim_x) ? 0 : (unsigned int) dim_x; //
  // actual dim to be operated on

  /* Step 1: do conditional in current dim */
  base::Grid* g_out = nullptr;
  base::DataVector* a_out = new base::DataVector(1);
  op_factory::createOperationDensityConditional(*g_in)->doConditional(
      *a_in, g_out, *a_out, static_cast<unsigned int>(op_dim), coords1d->get(curr_dim));

  // move on to next dim
  curr_dim = (curr_dim + 1) % dims;
  op_dim = (op_dim + 1) % g_out->getDimension();

  /* Step 2: draw a sample in next dim */
  double y = 0;

  if (g_out->getDimension() > 1) {
    // Marginalize to next dimension
    base::Grid* g1d = nullptr;
    base::DataVector* a1d = nullptr;
    op_factory::createOperationDensityMargTo1D(*g_out)->margToDimX(a_out, g1d, a1d, op_dim);

    // Draw a sample in next dimension
    y = doTransformation1D(g1d, a1d, coords1d->get(curr_dim));
    delete g1d;
    delete a1d;

  } else {
    // skip Marginalize, directly draw a sample in next dimension
    y = doTransformation1D(g_out, a_out, coords1d->get(curr_dim));
  }

  /* Step 3: copy sample to output */
  cdfs1d->set(curr_dim, y);

  /* Step 4: sample in next dimension */
  if (g_out->getDimension() > 1)
    doTransformation_in_next_dim(g_out, a_out, op_dim, coords1d, cdfs1d, curr_dim);

  delete g_out;
  delete a_out;

  return;
}

double OperationRosenblattTransformationModPolyClenshawCurtis::
    doTransformation1D(base::Grid* grid1d, base::DataVector* alpha1d, double coord1d) {
  OperationRosenblattTransformation1DModPolyClenshawCurtis* opRosenblatt
    = static_cast<OperationRosenblattTransformation1DModPolyClenshawCurtis*>
    (op_factory::createOperationRosenblattTransformation1D(*grid1d));
  double y = opRosenblatt->doTransformation1D(alpha1d, coord1d);
  delete opRosenblatt;
  if (y == 0) {
    std::cout << "Rosenblatt y=0" << std::endl;
    std::cout << alpha1d->toString() << std::endl;
    std::cout << "coord1d:" << coord1d << std::endl;
  }
  return y;
}  // end of compute_1D_cdf()
}  // namespace datadriven
}  // namespace sgpp
