// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/optimization/TrisectionOptimizer.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/globaldef.hpp>

#include <cmath>
#include <iostream>
#include <vector>

const double tolerance = 1e-8;

BOOST_AUTO_TEST_CASE(testTrisectionOptimizer) {
  double minValue = 0.123456789;
  sgpp::combigrid::SingleFunction f(
      [minValue](double x) { return (x - minValue) * (x - minValue); });
  sgpp::combigrid::TrisectionOptimizer trOpt(f);

  sgpp::combigrid::Stopwatch watch;

  auto result = trOpt.refine(sgpp::combigrid::OptimizationGuess::initial(0.0, 1.0, f));

  watch.log();
  std::cout << "Difference: " << fabs(result.b - minValue) << "\n";

  BOOST_CHECK_CLOSE(result.b, minValue, tolerance);
}
