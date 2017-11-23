// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iostream>

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridGridBasedEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>
#include <sgpp/globaldef.hpp>

using sgpp::combigrid::MultiIndex;

double atanModel(sgpp::base::DataVector const& v) {
  //  return std::atan(50 * (v[0] - .35));
  //  return v[0] * v[1] * v[1];
  return std::atan(50 * (v[0] - .35)) + M_PI / 2 + 4 * std::pow(v[1], 3) +
         std::exp(v[0] * v[1] - 1);
}

// This data was created by SGpp/combigrid/tests/createVarianceData.py
// It represents the variance calculated for each level with levelsum <= 7 for the function
// arctan(50.0 * (x[0] - .35)) + pi / 2.0 + 4.0 * x[1] ** 3 + exp(x[0] * x[1] - 1.0)
struct VarianceTestData {
  std::vector<MultiIndex> levels{
      MultiIndex{0, 0}, MultiIndex{0, 1}, MultiIndex{0, 2}, MultiIndex{0, 3}, MultiIndex{0, 4},
      MultiIndex{0, 5}, MultiIndex{0, 6}, MultiIndex{0, 7}, MultiIndex{1, 0}, MultiIndex{1, 1},
      MultiIndex{1, 2}, MultiIndex{1, 3}, MultiIndex{1, 4}, MultiIndex{1, 5}, MultiIndex{1, 6},
      MultiIndex{2, 0}, MultiIndex{2, 1}, MultiIndex{2, 2}, MultiIndex{2, 3}, MultiIndex{2, 4},
      MultiIndex{2, 5}, MultiIndex{3, 0}, MultiIndex{3, 1}, MultiIndex{3, 2}, MultiIndex{3, 3},
      MultiIndex{3, 4}, MultiIndex{4, 0}, MultiIndex{4, 1}, MultiIndex{4, 2}, MultiIndex{4, 3},
      MultiIndex{5, 0}, MultiIndex{5, 1}, MultiIndex{5, 2}, MultiIndex{6, 0}, MultiIndex{6, 1},
      MultiIndex{7, 0}};
  std::vector<double> variances{
      0,        1.701157, 1.437021, 1.437021, 1.437021, 1.437021, 1.437021, 1.437021, 1.080108,
      2.816791, 2.549883, 2.549881, 2.549880, 2.549880, 2.549880, 2.545726, 4.286104, 4.019262,
      4.019260, 4.019259, 4.019259, 1.858681, 3.597641, 3.330804, 3.330801, 3.330801, 1.985492,
      3.725129, 3.458287, 3.458285, 1.979975, 3.719659, 3.452817, 1.980263, 3.719941, 1.980269};
};

double BSplineVariance(sgpp::combigrid::MultiIndex level) {
  size_t numDimensions = 2;
  size_t degree = 3;
  sgpp::combigrid::MultiFunction func(atanModel);
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());

  sgpp::combigrid::EvaluatorConfiguration evalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineScalarProduct, degree);

  sgpp::combigrid::CombiEvaluators::MultiCollection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(evalConfig));

  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);
  bool exploitNesting = false;
  auto summationStrategyType = sgpp::combigrid::FullGridSummationStrategyType::VARIANCE;

  auto storage = std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage>(
      new sgpp::combigrid::CombigridTreeStorage(pointHierarchies, exploitNesting));

  std::shared_ptr<sgpp::combigrid::AbstractFullGridEvaluator<sgpp::combigrid::FloatArrayVector>>
      fullGridEval = std::make_shared<
          sgpp::combigrid::FullGridGridBasedEvaluator<sgpp::combigrid::FloatArrayVector>>(
          storage, evaluators, pointHierarchies, gf, summationStrategyType);

  auto result = fullGridEval->eval(level);
  double res = result[0].value();
  return res;
}

BOOST_AUTO_TEST_SUITE(testBsplines)

/*
 * This test calculates the variance for the test function atanModel (see above) and compares it to
 * precalculated data.
 *
 * ToDo (rehmemk) why is there an offset of 1e-7 between python data and quadrature data. Shouldn't
 * it be much smaller?
 */

BOOST_AUTO_TEST_CASE(testVarianceOnLevel) {
  std::cout << "Testing B spline variance calculation  subgridwise on single levels" << std::endl;
  VarianceTestData varianceTestData;
  double tolerance = 1e-6;
  double varianceError;
  for (size_t i = 0; i < varianceTestData.levels.size(); i++) {
    MultiIndex level = varianceTestData.levels[i];
    double bSplineVariance = BSplineVariance(level);
    varianceError = std::fabs(bSplineVariance - varianceTestData.variances[i]);
    std::cout << "level: " << level[0] << " " << level[1] << "|  error:  " << varianceError
              << std::endl;
    BOOST_CHECK_SMALL(varianceError, tolerance);
  }
}

BOOST_AUTO_TEST_SUITE_END()
