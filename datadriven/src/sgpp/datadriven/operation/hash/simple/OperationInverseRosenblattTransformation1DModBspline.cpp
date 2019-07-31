// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/tools/HermiteBasis.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DModBspline.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DModBspline.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp_datadriven.hpp>
#include <sgpp_optimization.hpp>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {
/**
 * WARNING: the grid must be a 1D grid!
 */
OperationInverseRosenblattTransformation1DModBspline::
    OperationInverseRosenblattTransformation1DModBspline(base::Grid* grid)
    : sum(0.0), quadOrder(0), grid(grid) {}

OperationInverseRosenblattTransformation1DModBspline::
    ~OperationInverseRosenblattTransformation1DModBspline() {}

void OperationInverseRosenblattTransformation1DModBspline::init(base::DataVector* alpha1d) {
  patch_areas.clear();
  is_negative_patch.clear();
  ordered_grid_points.clear();
  patch_functions.clear();
  coord_cdf.clear();
  sum = 0;
  opEval = std::unique_ptr<base::OperationEval>(op_factory::createOperationEvalNaive(*grid));

  base::DataVector coord(1);
  std::multimap<double, double>::iterator it1;
  base::GridStorage* gs = &this->grid->getStorage();
  double area = 0.0;
  double right_coord = 0.0, right_function_value = 0.0;
  size_t p = dynamic_cast<sgpp::base::ModBsplineGrid*>(grid)->getDegree();
  quadOrder = (p + 1) / 2;

  gauss.getLevelPointsAndWeightsNormalized(quadOrder, gauss_coordinates, weights);

  // need an ordered list of the grid points
  for (size_t i = 0; i < gs->getSize(); i++) {
    coord[0] = gs->getPoint(i).getStandardCoordinate(0);
    ordered_grid_points.push_back(coord[0]);
    // inserting dummy points into coord_cdf to make it have the right size
    coord_cdf.insert(std::pair<double, double>(coord[0], 0.0));
  }
  ordered_grid_points.push_back(0.0);
  ordered_grid_points.push_back(1.0);
  std::sort(ordered_grid_points.begin(), ordered_grid_points.end());

  coord_cdf.insert(std::pair<double, double>(0.0, 0.0));
  coord_cdf.insert(std::pair<double, double>(1.0, 1.0));

  double left_coord = 0.0;
  coord[0] = 0.0;
  double left_function_value = opEval->eval(*alpha1d, coord);
  bool negative_start = left_function_value < 0;
  left_function_value = std::max(0.0, left_function_value);
  for (size_t i = 1; i < ordered_grid_points.size(); i++) {
    coord[0] = ordered_grid_points[i];
    double eval_res = opEval->eval(*alpha1d, coord);

    double gaussQuadSum = 0.;
    double left = left_coord;
    double scaling = coord[0] - left;
    // this will always be initialized to false except if we are in the first patch
    // and the first value equals 0
    bool negative_value_encountered = (i == 1 && negative_start);
    for (size_t c = 0; c < quadOrder; c++) {
      coord[0] = left + scaling * gauss_coordinates[c];
      double value = opEval->eval(*alpha1d, coord);
      gaussQuadSum += weights[c] * value;
      negative_value_encountered = negative_value_encountered || (value < 0);
    }
    area = gaussQuadSum * scaling;

    if (negative_value_encountered || eval_res < 0) {
      // make sure that the cdf is monotonically increasing
      // WARNING: THIS IS A HACK THAT OVERCOMES THE PROBLEM
      // OF NON POSITIVE DENSITY

      // we look for the next grid point with pdf(x) >= 0
      size_t j;
      for (j = i; j < ordered_grid_points.size() - 1; j++) {
        right_coord = ordered_grid_points[j];
        coord[0] = right_coord;
        right_function_value = opEval->eval(*alpha1d, coord);
        if (right_function_value >= 0 && right_function_value != left_function_value) break;
      }

      right_coord = ordered_grid_points[j];
      if (j == ordered_grid_points.size() - 1) {
        if (left == 0)
          right_function_value = 1;
        else
          right_function_value = 0;
      }

      // get last function value and coordinate with pdf(x) >= 0
      // perform montonic cubic interpolation based on:
      // https://en.wikipedia.org/wiki/Monotone_cubic_interpolation

      // cubic interpolation -> use three points:
      // - the last with pdf(x) >= 0 (left from negative interval)
      // - the first with pdf(x) >= 0 (right from negative inteval)
      // - the right neighbor of the right one
      double secants[2];
      double tangents[3];
      double function_values[3];
      function_values[0] = left_function_value;
      function_values[1] = right_function_value;
      secants[0] = (right_function_value - left_function_value) / (right_coord - left_coord);
      tangents[0] = secants[0];
      if (j != ordered_grid_points.size() - 1) {
        coord[0] = ordered_grid_points[j + 1];
        function_values[2] = opEval->eval(*alpha1d, coord);
      } else {
        // if j is the last grid point choose the next one with the same step size
        // and set it's function value to zero
        // and set it's function value to the last value
        coord[0] = 1 + ordered_grid_points[j] - ordered_grid_points[j - 1];
        function_values[2] = right_function_value;
      }
      secants[1] = (function_values[2] - function_values[1]) / (coord[0] - ordered_grid_points[j]);
      tangents[2] = secants[1];
      // secants left and right of current point
      if (secants[0] == 0 || secants[1] == 0)
        // if one of the secants is zero
        tangents[1] = 0;
      else if ((secants[0] > 0 && secants[1] < 0) || (secants[0] < 0 && secants[1] > 0))
        // if the secants dont have the same sign
        tangents[1] = 0;
      else
        tangents[1] = (secants[0] + secants[1]) / 2;

      // correction to make the interpolation strict monotonic
      for (size_t c = 0; c < 2; c++) {
        double alpha = tangents[c] / secants[c];
        double beta = tangents[c + 1] / secants[c];
        if (alpha * alpha + beta * beta > 9) {
          double tau = 3. / std::sqrt(alpha * alpha + beta * beta);
          tangents[c] = tau * alpha * secants[c];
          tangents[c + 1] = tau * beta * secants[c];
        }
      }
      // interpolation that can be evaluated between left_coord and right_coord

      std::function<double(double)> interpolation = [right_coord, left_coord, left_function_value,
                                                     right_function_value,
                                                     tangents](double x) -> double {
        double h = right_coord - left_coord;
        double t = (x - left_coord) / h;
        return left_function_value * base::HermiteBasis::h_0_0(t) +
               h * tangents[0] * base::HermiteBasis::h_1_0(t) +
               +right_function_value * base::HermiteBasis::h_0_1(t) +
               +h * tangents[1] * base::HermiteBasis::h_1_1(t);
      };

      for (; i <= j; i++) {
        coord[0] = ordered_grid_points[i];
        double gaussQuadSum = 0.;
        double left = left_coord;
        double scaling = coord[0] - left;
        for (size_t c = 0; c < quadOrder; c++) {
          coord[0] = left + scaling * gauss_coordinates[c];
          gaussQuadSum += weights[c] * interpolation(coord[0]);
        }

        area = gaussQuadSum * scaling;
        sum += area;
        patch_areas.push_back(area);
        is_negative_patch.push_back(true);
        patch_functions.push_back(interpolation);
        left_coord = ordered_grid_points[i];
      }
      --i;
      left_function_value = right_function_value;
    } else {
      // use the (positive) result of the Gauss-Quadrature
      sum += area;
      left_coord = left + scaling;
      left_function_value = eval_res;
      patch_areas.push_back(area);
      is_negative_patch.push_back(false);
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
}

double OperationInverseRosenblattTransformation1DModBspline::sample(base::DataVector* alpha1d,
                                                                    double coord1d) {
  if (coord1d == 0.0) return 0.0;
  if (sum == 0) return 0;
  base::DataVector coord(1);
  std::multimap<double, double>::iterator it1;
  size_t patch_nr = 0;
  int negative_patch_counter = 0;
  for (it1 = coord_cdf.begin(); it1 != coord_cdf.end(); ++it1) {
    if (it1->first == coord1d)
      return it1->second;
    else if (it1->first >= coord1d)
      break;
    if (is_negative_patch[patch_nr]) {
      // std::cout << "patch " << patch_nr << "is negative" << std::endl;
      ++negative_patch_counter;
    }
    ++patch_nr;
  }
  --patch_nr;
  --negative_patch_counter;
  // std::cout << "negative_patch_counter: " << negative_patch_counter << std::endl;
  // std::cout << "patch_nr: " << patch_nr << std::endl;
  // std::cout << "PAs size: " << patch_areas.size() << std::endl;
  // std::cout << "PFs size : " << patch_functions.size() << std::endl;
  --it1;
  double gaussQuadSum = 0.;
  double left = it1->first;
  double scaling = coord1d - left;
  if (is_negative_patch[patch_nr]) {
    for (size_t c = 0; c < quadOrder; c++) {
      coord[0] = left + scaling * gauss_coordinates[c];
      gaussQuadSum += weights[c] * patch_functions[negative_patch_counter](coord[0]);
    }
  } else {
    for (size_t c = 0; c < quadOrder; c++) {
      coord[0] = left + scaling * gauss_coordinates[c];
      gaussQuadSum += weights[c] * opEval->eval(*alpha1d, coord);
    }
  }
  return it1->second + (gaussQuadSum * scaling) / sum;
}

double OperationInverseRosenblattTransformation1DModBspline::doTransformation1D(
    base::DataVector* alpha1d, double coord1d) {
  init(alpha1d);
  // std::cout << "PFs size after exit: " << patch_functions.size() << std::endl;
  std::function<double(const base::DataVector&)> optFunc = [this, coord1d, alpha1d](
      const base::DataVector& x) -> double {
    double F_x = sample(alpha1d, x[0]);
    return (F_x - coord1d) * (F_x - coord1d);
  };
  base::Printer::getInstance().disableStatusPrinting();
  base::WrapperScalarFunction f(1, optFunc);
  optimization::optimizer::NelderMead nelderMead(f);
  nelderMead.optimize();
  return nelderMead.getOptimalPoint()[0];
}
}  // namespace datadriven
}  // namespace sgpp
