// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DPoly.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>

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
OperationRosenblattTransformation1DPoly::OperationRosenblattTransformation1DPoly(
    base::Grid* grid)
    : grid(grid) {}

OperationRosenblattTransformation1DPoly::~OperationRosenblattTransformation1DPoly() {}

double OperationRosenblattTransformation1DPoly::doTransformation1D(base::DataVector* alpha1d,
                                                                  double coord1d) {
  if (coord1d == 0.0)
    return 0.0;
  /***************** STEP 1. Compute CDF  ********************/

  // compute PDF, sort by coordinates
  std::multimap<double, double> coord_cdf;
  std::vector<double> ordered_grid_points;
  std::multimap<double, double>::iterator it1;

  base::GridStorage* gs = &this->grid->getStorage();
  std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*grid));
  base::DataVector coord(1);
  base::DataVector gauss_coordinates;
  base::DataVector weights;
  base::GaussLegendreQuadRule1D gauss;
  std::vector<double> patch_areas;
  size_t p = dynamic_cast<sgpp::base::PolyGrid*>(grid)->getDegree();
  const size_t quadOrder =  (p + 1) / 2;
  gauss.getLevelPointsAndWeightsNormalized(quadOrder, gauss_coordinates, weights);

  double right_coord, right_function_value;
  std::cout << "size:" << gs->getSize() << std::endl;
  double sum = 0.0 , area = 0.0;

  // need an ordered list of the grid points
  for (size_t i = 0; i < gs->getSize(); i++) {
    coord[0] = gs->getPoint(i).getStandardCoordinate(0);
    ordered_grid_points.push_back(coord[0]);
    // inserting dummy points into coord_cdf to make it have the right size
    coord_cdf.insert(std::pair<double, double>(coord[0], 0.0));
  }

  ordered_grid_points.push_back(0.0);
  ordered_grid_points.push_back(1.0);

  coord_cdf.insert(std::pair<double, double>(0.0, 0.0));
  coord_cdf.insert(std::pair<double, double>(1.0, 1.0));
  std::sort(ordered_grid_points.begin(), ordered_grid_points.end());

  double left_coord = 0.0;
  double left_function_value = 0.0;
  for (size_t i = 1; i < ordered_grid_points.size(); i++) {
    coord[0] = ordered_grid_points[i];
    double eval_res = opEval->eval(*alpha1d, coord);
    std::cout << "pdf("
             << coord[0] << ")=" << eval_res << std::endl;
    double gaussQuadSum = 0.;
    double left = left_coord;
    double scaling = coord[0] - left;
    for (size_t c = 0; c < quadOrder; c++) {
      coord[0] = left + scaling * gauss_coordinates[c];
      gaussQuadSum += weights[c] * opEval->eval(*alpha1d, coord);
    }
    area = gaussQuadSum * scaling;
    std::cout << "from " << left_coord << " to " << left+scaling << std::endl;
    std::cout << "area:" << area << std::endl;

    if (area < 0) {
      // make sure that the cdf is monotonically increasing
      // WARNING: THIS IS A HACK THAT OVERCOMES THE PROBLEM
      // OF NON POSITIVE DENSITY
      std::cerr << "warning: negative pdf value encountered pdf("
               << coord[0] << ")=" << eval_res
               << std::endl;

      // we look for the next grid point with pdf(x) >= 0
      size_t j;
      for (j = i; j < ordered_grid_points.size(); j++) {
        std::cout << "j:" << j << std::endl;
        right_coord = ordered_grid_points[j];
        coord[0] = right_coord;
        right_function_value = opEval->eval(*alpha1d, coord);
        if (right_function_value >= 0)
          break;
      }
      std::cout << "Found j: " << j << std::endl;
      std::cout << right_coord << ";" << right_function_value << std::endl;
      // get last function value and coordinate with pdf(x) > 0
      // interpolate every value in between linearly
      for (; i < j; i++) {
        coord[0] = ordered_grid_points[i];
        std::cout << "interpolating i:" << i << std::endl;
        eval_res = left_function_value
          + (right_function_value + left_function_value) / (right_coord - left_coord)
          * (coord[0] - left_coord);
        std::cout << "For x=" << coord[0] << "interp: " << eval_res << std::endl;
        area = (eval_res + left_function_value) / 2 * (coord[0] - left_coord);
        sum += area;
        std::cout << "from " << left_coord << " to " << coord[0] << std::endl;
        std::cout << "area:" << area << std::endl;
        patch_areas.push_back(area);
        left_function_value = eval_res;
        left_coord = coord[0];
      }
      area = (right_function_value + left_function_value) / 2 * (right_coord - left_coord);
      sum += area;
      std::cout << "from " << left_coord << " to " << right_coord << std::endl;
      std::cout << "area:" << area << std::endl;
      patch_areas.push_back(area);
      left_coord = right_coord;
      left_function_value = right_function_value;
    } else {
      // use the (positive) result of the Gauss-Quadrature
      sum += area;
      left_coord = left + scaling;
      left_function_value = eval_res;
      patch_areas.push_back(area);
    }
  }

  // compute CDF
  double tmp_sum;
  unsigned int i = 0;

  for (it1 = coord_cdf.begin(); it1 != coord_cdf.end(); ++it1) {
    tmp_sum = 0.0;
    for (unsigned int j = 0; j < i; ++j) tmp_sum += patch_areas[j];
    ++i;
    it1->second = tmp_sum / sum;
  }

  /***************** STEP 1. Done  ********************/

  /***************** STEP 2. Sampling  ********************/
  double y, x1, x2, y1, y2;

  std::cout << "Areas: " << std::endl;
  for (size_t i = 0; i < patch_areas.size(); i++) {
    std::cout << patch_areas[i] << std::endl;
  }
  std::cout << "Size areas: " << patch_areas.size() << std::endl;
  std::cout << "Size cdf: " << coord_cdf.size() << std::endl;
  std::cout << "coord cdf: " << std::endl;
  for (it1 = coord_cdf.begin(); it1 != coord_cdf.end(); ++it1) {
    std::cout << it1->second << std::endl;
  }
  // find cdf interval
  for (it1 = coord_cdf.begin(); it1 != coord_cdf.end(); ++it1) {
    if (it1->first >= coord1d) break;
  }

  x2 = it1->first;
  y2 = it1->second;
  --it1;
  x1 = it1->first;
  y1 = it1->second;
  // find x (linear interpolation): (y-y1)/(x-x1) = (y2-y1)/(x2-x1)
  y = (y2 - y1) / (x2 - x1) * (coord1d - x1) + y1;

  /***************** STEP 2. Done  ********************/
  return y;
}  // end of compute_1D_cdf()
}  // namespace datadriven
}  // namespace sgpp
