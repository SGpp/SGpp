// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DLinear.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>
#include <map>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
#include <algorithm>

namespace sgpp {
namespace datadriven {
/**
 * WARNING: the grid must be a 1D grid!
 */
OperationRosenblattTransformation1DLinear::OperationRosenblattTransformation1DLinear(
    base::Grid* grid)
    : grid(grid) {}

OperationRosenblattTransformation1DLinear::~OperationRosenblattTransformation1DLinear() {}

double OperationRosenblattTransformation1DLinear::doTransformation1D(base::DataVector* alpha1d,
                                                                     double coord1d) {
  /***************** STEP 1. Compute CDF  ********************/

  // compute PDF, sort by coordinates
  std::multimap<double, double> coord_pdf, coord_cdf;
  std::multimap<double, double>::iterator it1, it2;

  base::GridStorage* gs = &this->grid->getStorage();
  std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*(this->grid)));
  base::DataVector coord(1);

  for (size_t i = 0; i < gs->getSize(); i++) {
    coord[0] = gs->getPoint(i).getStandardCoordinate(0);
    coord_pdf.insert(std::pair<double, double>(coord[0], opEval->eval(*alpha1d, coord)));
    coord_cdf.insert(std::pair<double, double>(coord[0], i));
  }

  // include values at the boundary [0,1]
  coord_pdf.insert(std::pair<double, double>(0.0, 0.0));
  coord_pdf.insert(std::pair<double, double>(1.0, 0.0));
  coord_cdf.insert(std::pair<double, double>(0.0, 0.0));
  coord_cdf.insert(std::pair<double, double>(1.0, 1.0));

  // Composite rule: trapezoidal (b-a)/2 * (f(a)+f(b))
  it1 = coord_pdf.begin();
  it2 = coord_pdf.begin();
  std::vector<double> tmp;
  tmp.push_back(0.0);
  double sum = 0.0, area = 0.0;

  for (++it2; it2 != coord_pdf.end(); ++it2) {
    // (*it).first : the coordinate
    // (*it).second : the function value
    area = (it2->first - it1->first) * (it1->second + it2->second) / 2;

    // make sure that the cdf is monotonically increasing
    // WARNING: THIS IS A HACK THAT OVERCOMES THE PROBLEM
    // OF NON POSITIVE DENSITY
    if (area < 0) {
      std::cerr << "warning: negative area encountered " << (*it1).second << ", " << (*it2).second
                << std::endl;
      area = 0;
    }

    sum += area;
    tmp.push_back(area);
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
