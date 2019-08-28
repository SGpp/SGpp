// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/optimization/MixedOptimizer.hpp>
#include <sgpp/combigrid/optimization/NewtonOptimizer.hpp>
#include <sgpp/combigrid/optimization/TrisectionOptimizer.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/globaldef.hpp>

#include <cmath>
#include <iostream>
#include <vector>

const double tolerance = 1e-5;

BOOST_AUTO_TEST_CASE(testTrisectionOptimizer) {
  double minValue = 0.123456789;
  // The +1 makes optimization a lot more imprecise because of cancelling effects around the minimum
  sgpp::combigrid::SingleFunction f(
      [minValue](double x) { return (x - minValue) * (x - minValue) + 1; });
  sgpp::combigrid::TrisectionOptimizer trOpt(f);

  sgpp::combigrid::Stopwatch watch;

  auto result = trOpt.minimize(sgpp::combigrid::OptimizationGuess::initial(0.0, 1.0, f), 50, 10);

  // watch.log();
  // std::cout << "Difference: " << fabs(result.b - minValue) << "\n";

  BOOST_CHECK_CLOSE(result.b, minValue, tolerance);
}

BOOST_AUTO_TEST_CASE(testNewtonOptimizer) {
  double minValue = 0.123456789;
  sgpp::combigrid::SingleFunction f(
      [minValue](double x) { return (x - minValue) * (x - minValue) + 1; });

  sgpp::combigrid::NewtonOptimizer nOpt(f);

  sgpp::combigrid::Stopwatch watch;

  auto result = nOpt.minimize(sgpp::combigrid::OptimizationGuess::initial(0.0, 1.0, f), 3);

  // watch.log();
  // std::cout << "Difference: " << fabs(result.b - minValue) << "\n";

  BOOST_CHECK_CLOSE(result.b, minValue, tolerance);
}

BOOST_AUTO_TEST_CASE(testMixedOptimizer) {
  double minValue = M_PI / 2.0;
  sgpp::combigrid::SingleFunction f([](double x) { return -sin(x); });

  sgpp::combigrid::MixedOptimizer opt(f);

  sgpp::combigrid::Stopwatch watch;

  auto result = opt.minimize(sgpp::combigrid::OptimizationGuess::initial(0.0, 2.0, f));

  // watch.log();
  // std::cout << "Difference: " << fabs(result.b - minValue) << "\n";

  BOOST_CHECK_CLOSE(result.b, minValue, tolerance);
}
