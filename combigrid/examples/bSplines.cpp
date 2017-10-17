// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/quadrature/sampling/NaiveSampleGenerator.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

double f(sgpp::base::DataVector const& v) {
  //  return v[0];  // * v[0] * v[0];
  //  return sin(1. / (1. + v[0] * v[0])) * v[1];
  //  return v[0];
  return std::atan(50 * (v[0] - .35)) + M_PI / 2 + 4 * std::pow(v[1], 3) +
         std::exp(v[0] * v[1] - 1);
}

/*
 * Borehole function models water flow through a borehole. Response is water flow rate in m^3/year
 * https://www.sfu.ca/~ssurjano/borehole.html
 *
 * Input:
 * rw in [0.05,0.15] radius of borehole (m)
 * r in [100,50000] radius of influence (m)
 * Tu in [63070,115600] transmissivity of upper aquifier (m^2/year)
 * Hu in [990,1110] potentiometric head of upper aquifier (m)
 * Tl in [63.1 , 116] transmissivity of lower aquifier(m^2/year)
 * Hl in [700,820] potentiometric head of lower aquifier (m)
 * L in [1120, 1680] length of borehole (m)
 * Kw in [9855, 12045] hydraulic conductivity of borehole (m/year)
 */

double borehole(double rw, double r, double Tu, double Hu, double Tl, double Hl, double L,
                double Kw) {
  double dominator = 2 * M_PI * Tu * (Hu - Hl);
  double logterm = log(r / rw);
  double denominator = 1 + (2 * L * Tu / (logterm * rw * rw * Kw)) + Tu / Tl;
  return dominator / (logterm * denominator);
}

void interpolate(size_t maxlevel, double& max_err, double& L2_err) {
  size_t numDimensions = 2;
  size_t degree = 5;
  sgpp::combigrid::MultiFunction func(f);
  auto operation =
      sgpp::combigrid::CombigridOperation::createExpUniformBoundaryBsplineInterpolation(
          numDimensions, func, degree);

  //  auto operation =
  //      sgpp::combigrid::CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(
  //          numDimensions, func);

  //  auto operation =
  //  sgpp::combigrid::CombigridOperation::createExpUniformBoundaryLinearInterpolation(
  //      numDimensions, func);

  double diff = 0.0;
  // generator generates num_points random points in [0,1]^numDimensions
  size_t num_points = 1000;
  sgpp::quadrature::NaiveSampleGenerator generator(numDimensions);
  sgpp::base::DataVector p(numDimensions, 0);

  for (size_t i = 0; i < num_points; i++) {
    generator.getSample(p);
    diff = fabs(operation->evaluate(maxlevel, p) - f(p));
    max_err = (diff > max_err) ? diff : max_err;
    L2_err += pow(diff, 2);
  }
  L2_err = sqrt(L2_err / static_cast<double>(num_points));

  std::cout << "# grid points: " << operation->numGridPoints() << " ";
}

int main() {
  sgpp::base::SGppStopwatch watch;
  watch.start();
  size_t minLevel = 0;
  size_t maxLevel = 8;
  std::vector<double> maxErr(maxLevel + 1, 0);
  std::vector<double> L2Err(maxLevel + 1, 0);
  for (size_t l = minLevel; l < maxLevel + 1; l++) {
    interpolate(l, maxErr[l], L2Err[l]);
    std::cout << "level: " << l << " max err " << maxErr[l] << " L2 err " << L2Err[l] << std::endl;
  }
  std::cout << " Total Runtime: " << watch.stop() << " s" << std::endl;
  return 0;
}
