// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "bSplines.hpp"

int main() {
  //  size_t numDimensions = 1;
  size_t degree = 3;
  //  size_t level = 2;
  sgpp::combigrid::MultiIndex oneLevel = {2, 2};

  // Interpolation
  //  sgpp::base::SGppStopwatch watch;
  //  watch.start();
  //  size_t minLevel = 2;
  //  size_t maxLevel = 2;
  //
  //  std::vector<double> maxErr(maxLevel + 1, 0);
  //  std::vector<double> L2Err(maxLevel + 1, 0);
  //  for (size_t l = minLevel; l < maxLevel + 1; l++) {
  //    interpolate(l, numDimensions, degree, maxErr[l], L2Err[l]);
  //    //    std::cout << "level: " << l << " max err " << maxErr[l] << " L2 err " << L2Err[l] <<
  //    //    std::endl;
  //  }

  //  std::cout << " Total Runtime: " << watch.stop() << " s" << std::endl;

  // Integration
  //  double integral = integrate(level, numDimensions, degree);
  //  std::cout << "integral:  " << integral << std::endl;

  // Integrate basis functions
  //  std::vector<double> integrals = integrateBasisFunctions(level, numDimensions, degree);
  //  std::cout << "------------------------------------" << std::endl;
  //  for (size_t i = 0; i < integrals.size(); i++) {
  //    std::cout << integrals[i] << " ";
  //  }

  //  Interpolate in one direction and integrate in the other
  //  double res = interpolate_and_integrate(level, numDimensions, degree);
  //  std::cout << res << std::endl;

  // Calculate integral of func^2
  //  double isqu = integrateSquare(level, numDimensions, degree);
  //  std::cout << "integral f^2 : " << isqu << std::endl;

  // Calculate variances on subgrids
  //  double var = variance(level, numDimensions, degree);
  //  std::cout << var << std::endl;

  //
  //  double value = interpolateOneLevel(oneLevel, numDimensions, degree);
  //  std::cout << "result: " << value << std::endl;

  for (size_t auxLevel = 0; auxLevel < 9; auxLevel++) {
    oneLevel = {auxLevel};
    double value = integrateOneLevel(oneLevel, degree);
    std::cout << auxLevel << " | integral: " << value << std::endl;
  }

  return 0;
}
