// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/*
 * This example illustrates how to use the PolynomialChaosExpansion Class.
 */
#include <cmath>
#include <iomanip>
#include <sgpp/datadriven/tools/PolynomialChaosExpansion.hpp>
#include <sgpp_base.hpp>
#include <string>
#include <utility>
#include <vector>

// define target function e
double e(const sgpp::base::DataVector& vec) {
  return 1 - ((4 * vec[1]) / (5 * 225 * vec[0])) -
         ((vec[2] * vec[2]) / (25 * 225 * vec[0] * vec[0]));
}

int main() {
  std::cout << std::fixed;
  std::cout << std::setprecision(9);
  // construct a DistributionsVector and initialize it
  sgpp::base::DistributionsVector dists;
  auto dist1 = std::make_shared<sgpp::base::DistributionLogNormal>(5, .5);
  auto dist2 = std::make_shared<sgpp::base::DistributionNormal>(2000, 400);
  auto dist3 = std::make_shared<sgpp::base::DistributionNormal>(500, 100);
  dists.push_back(dist1);
  dists.push_back(dist2);
  dists.push_back(dist3);
  /*
   * construct the PCE object with the target function, the expansion order and the
   * DistributionsVector
   */
  sgpp::datadriven::PolynomialChaosExpansion pce =
      sgpp::datadriven::PolynomialChaosExpansion(e, 3, dists);
  /*
   * now the methods of the PCE object are shown
   *
   * the first parameter is the number of grid points that should be used
   * for regular sparse grids the smallest grid with at least this many grid points is used
   *
   * the second parameter in the following methods is a binary flag, signifying whether to use
   * adaptive sparse grids to calculate the coefficients
   *
   * Note that the coefficients are stored once computed, to compute them with different parameters
   * call pce.clearCoefficients()
   */
  pce.calculateCoefficients(800, true);
  std::cout << "-----------------------------------------------------------------------------------"
            << '\n';
  std::cout << "L2-Error: " << pce.getL2Error(800, true) << '\n';
  std::cout << "Expansion mean: " << pce.getMean(800, true) << '\n';
  std::cout << "Expansion variance: " << pce.getVariance(800, true) << '\n';

  std::cout << "Coefficients of the expansion: " << '\n';
  auto coeffs = pce.getCoefficients();
  for (auto entry : coeffs) {
    std::cout << entry << ", ";
  }
  std::cout << "-----------------------------------------------------------------------------------"
            << '\n';
}
