// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityConditional.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensitySampling1D.hpp>
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

void OperationRosenblattTransformationLinear::doTransformation(base::DataVector* alpha,
                                                               base::DataMatrix* points,
                                                               base::DataMatrix* pointscdf) {
  size_t dim_start = 0;
  size_t num_dims = this->grid->getDimension();
  size_t bucket_size = points->getNrows() / num_dims + 1;
  base::DataVector* coords1d = new base::DataVector(points->getNcols());
  base::DataVector* cdfs1d = new base::DataVector(points->getNcols());

  // 1. marginalize to dim_start
  base::Grid* g1d = NULL;
  base::DataVector* a1d = NULL;
  std::unique_ptr<OperationDensityMargTo1D> marg1d(
      op_factory::createOperationDensityMargTo1D(*this->grid));
  marg1d->margToDimX(alpha, g1d, a1d, dim_start);

  // 2. 1D transformation on dim_start
  double y = 0;

  // 3. for every sample do...
  // #pragma omp parallel
  //  {
  // #pragma omp for schedule(dynamic)
  for (size_t i = 0; i < points->getNrows(); i++) {
    //        OperationTransformation1D* trans1d =
    //        op_factory::createOperationRosenblattTransformation1D(*g1d);
    //        y = trans1d->doTransformation1D(a1d, points->get(i, dim_start));
    y = doTransformation1D(g1d, a1d, points->get(i, dim_start));
    pointscdf->set(i, dim_start, y);
    points->getRow(i, *coords1d);
    pointscdf->getRow(i, *cdfs1d);
    doTransformation_start_dimX(this->grid, alpha, dim_start, coords1d, cdfs1d);
    pointscdf->setRow(i, *cdfs1d);

    // change the starting dimension when the bucket_size is arrived
    // this distributes the error in the projection uniformly to all
    // dimensions and make it therefore stable
    if (((i + 1) % bucket_size) == 0 && (i + 1) < points->getNrows()) {
      dim_start++;
      // 1. marginalize to dim_start
      marg1d->margToDimX(alpha, g1d, a1d, dim_start);
    }
  }

  //  }

  delete g1d;
  delete a1d;
  delete coords1d;
  return;
}

void OperationRosenblattTransformationLinear::doTransformation(base::DataVector* alpha,
                                                               base::DataMatrix* points,
                                                               base::DataMatrix* pointscdf,
                                                               size_t dim_start) {
  // 1. marginalize to dim_start
  base::Grid* g1d = NULL;
  base::DataVector* a1d = NULL;
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
  return;
}

void OperationRosenblattTransformationLinear::doTransformation_start_dimX(
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

void OperationRosenblattTransformationLinear::doTransformation_in_next_dim(
    base::Grid* g_in, base::DataVector* a_in, size_t op_dim, base::DataVector* coords1d,
    base::DataVector* cdfs1d, size_t& curr_dim) {
  size_t dims = coords1d->getSize();  // total dimensions
  // unsigned int op_dim = curr_dim + 1;  // (curr_dim < dim_x) ? 0 : (unsigned int) dim_x; //
  // actual dim to be operated on

  /* Step 1: do conditional in current dim */
  base::Grid* g_out = NULL;
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
    base::Grid* g1d = NULL;
    base::DataVector* a1d = NULL;
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

double OperationRosenblattTransformationLinear::doTransformation1D(base::Grid* grid1d,
                                                                   base::DataVector* alpha1d,
                                                                   double coord1d) {
  /***************** STEP 1. Compute CDF  ********************/

  // compute PDF, sort by coordinates
  std::multimap<double, double> coord_pdf, coord_cdf;
  std::multimap<double, double>::iterator it1, it2, it3;

  base::GridStorage* gs = &grid1d->getStorage();
  std::unique_ptr<base::OperationEval> opEval = op_factory::createOperationEval(*(grid1d));
  base::DataVector coord(1);

  for (unsigned int i = 0; i < gs->getSize(); i++) {
    coord[0] = gs->get(i)->getCoord(0);
    coord_pdf.insert(std::pair<double, double>(coord[0], opEval->eval(*alpha1d, coord)));
    coord_cdf.insert(std::pair<double, double>(coord[0], i));
  }

  // include values at the boundary [0,1]
  coord_pdf.insert(std::pair<double, double>(0.0, 0.0));
  coord_pdf.insert(std::pair<double, double>(1.0, 0.0));
  coord_cdf.insert(std::pair<double, double>(0.0, 0.0));
  coord_cdf.insert(std::pair<double, double>(1.0, 1.0));

  // make sure that all the pdf values are positive
  // if not, interpolate between the closest positive neighbors
  it1 = coord_pdf.begin();
  it2 = coord_pdf.begin();

  it1->second = std::max(it1->second, 0.0);
  for (++it2; it2 != coord_pdf.end(); ++it2) {
    if (it2->second < 0.0) {
      // search for next right neighbor that has a positive function value
      it3 = it2;
      while (it3->second <= 0.0 && it3 != coord_pdf.end()) {
        it3++;
      }
      it2->second = (it1->second + it3->second) / 2.0;
    }

    it1 = it2;
  }

  // Composite rule: trapezoidal (b-a)/2 * (f(a)+f(b))
  it1 = coord_pdf.begin();
  it2 = coord_pdf.begin();
  std::vector<double> tmp;
  tmp.push_back(0.0);
  double sum = 0.0, area;

  for (++it2; it2 != coord_pdf.end(); ++it2) {
    // (*it).first : the coordinate
    // (*it).second : the function value
    area = ((*it2).first - (*it1).first) / 2 * ((*it1).second + (*it2).second);

    // make sure that the cdf is monotonically increasing
    // WARNING: THIS IS A HACK THAT OVERCOMES THE PROBLEM
    // OF NON POSITIVE DENSITY
    if (area < 0) {
      std::cerr << "warning: negative area encountered " << (*it1).second << ", " << (*it2).second
                << std::endl;
      area = 0;
    }

    tmp.push_back(area);
    sum += area;
    ++it1;
  }

  // compute CDF
  double tmp_sum;
  unsigned int i = 0;

  for (it1 = coord_cdf.begin(); it1 != coord_cdf.end(); ++it1) {
    tmp_sum = 0.0;

    for (unsigned int j = 0; j <= i; ++j) tmp_sum += tmp[j];

    ++i;
    (*it1).second = tmp_sum / sum;
  }

  tmp.clear();
  coord_pdf.clear();
  /***************** STEP 1. Done  ********************/

  /***************** STEP 2. Sampling  ********************/
  double y, x1, x2, y1, y2;

  // find cdf interval
  for (it1 = coord_cdf.begin(); it1 != coord_cdf.end(); ++it1) {
    if ((*it1).first >= coord1d) break;
  }

  x2 = (*it1).first;
  y2 = (*it1).second;
  --it1;
  x1 = (*it1).first;
  y1 = (*it1).second;
  // find x (linear interpolation): (y-y1)/(x-x1) = (y2-y1)/(x2-x1)
  y = (y2 - y1) / (x2 - x1) * (coord1d - x1) + y1;

  /***************** STEP 2. Done  ********************/
  return y;
}  // end of compute_1D_cdf()
}  // namespace datadriven
}  // namespace sgpp
