// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iostream>

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridGridBasedEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/quadrature/sampling/NaiveSampleGenerator.hpp>

using sgpp::combigrid::MultiIndex;

//--------------------- Preparations for testVarianceOneLevel---------------------------------------
double atanModel(sgpp::base::DataVector const& v) {
  //  return v[0] * v[0] * v[0] * v[1] * v[1] * v[1] * v[1];
  //  return std::atan(50 * (v[0] - .35));
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
      -0.000000000000003553, 1.701156850402689713, 1.437021035230856114, 1.437020728026304539,
      1.437020668119384226,  1.437020666114214862, 1.437020666048240969, 1.437020666045345507,
      1.080108355131402575,  2.816790782428100215, 2.549883454168121233, 2.549880760493220322,
      2.549880271385291053,  2.549880259715539665, 2.549880259378866754, 2.545726212447136483,
      4.286103665868155943,  4.019262025441856068, 4.019259816800595075, 4.019259376311369536,
      4.019259366026524560,  1.858681146494113534, 3.597641299962191397, 3.330803514248586339,
      3.330801307994530447,  3.330800874338359918, 1.985491657343512628, 3.725129287185531268,
      3.458287364652921525,  3.458285162235060994, 1.979975007603208326, 3.719658912837516596,
      3.452816643095617977,  1.980262657807472237, 3.719940998350237393, 1.980269030387413309};
};

struct VarianceDiagonalTestData {
  std::vector<MultiIndex> levels{MultiIndex{0, 0}, MultiIndex{1, 1}, MultiIndex{2, 2},
                                 MultiIndex{3, 3}, MultiIndex{4, 4}, MultiIndex{5, 5},
                                 MultiIndex{6, 6}, MultiIndex{7, 7}};
  double variance = 3.453102449971636290;
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

// --------------------------------Preparations for testMCVariance-- ---------------------------
namespace ishigami {
double tolerance = 5e-4;
size_t numDims = 3;
std::vector<double> pi_vec{M_PI, M_PI, M_PI};
sgpp::base::DataVector pi(pi_vec);
double pi_4 = M_PI * M_PI * M_PI * M_PI;
double a = 7.0;
double b = 0.1;
double variance = a * a / 8. + b * pi_4 / 5 + b * b * pi_4 * pi_4 / 18. + 0.5;
double mean = 3.5;
std::vector<double> sobolIndices{0.3138, 0.4424, 0.0, 0.0, 0.2436, 0.0, 0.0};
std::vector<double> totalSobolIndices{0.5574, 0.4424, 0.2436};
std::vector<double> bounds{0, 1};

double eval(sgpp::base::DataVector const& v) {
  // transform [0, 1] -> [-pi, pi]
  sgpp::base::DataVector x(v);
  x.mult(2 * M_PI);
  x.sub(ishigami::pi);

  // evaluate the Ishigami function
  return std::sin(x[0]) + ishigami::a * std::sin(x[1]) * std::sin(x[1]) +
         ishigami::b * x[2] * x[2] * x[2] * x[2] * std::sin(x[0]);
}
}  // namespace ishigami

namespace parabola {
double tolerance = 1e-13;
size_t numDims = 2;

double alpha1 = 5.0;
double beta1 = 4.0;

double alpha2 = 3.0;
double beta2 = 2.0;

double c1 = std::tgamma(alpha1 + beta1) / (std::tgamma(alpha1) * std::tgamma(beta1));
double c2 = std::tgamma(alpha2 + beta2) / (std::tgamma(alpha2) * std::tgamma(beta2));

double mean = c1 * c2 / 4725.0;
double variance = c1 * c1 * c1 * c2 * c2 * c2 / 75014100000.0 - 2 * c1 * c1 * c2 * c2 / 22325625.0 +
                  4.0 * c1 * c2 / 24255.0;
std::vector<double> bounds{0, 1};

double eval(sgpp::base::DataVector const& v) {
  double ans = 1.0;
  for (size_t idim = 0; idim < v.getSize(); idim++) {
    ans *= 4.0 * v[idim] * (1.0 - v[idim]);
  }
  return ans;
}
}  // namespace parabola

namespace debugfct {
size_t numDims = 1;
std::vector<double> bounds{0, 1};
double mean = 0.5;
double variance = 0.083333333333333333333;
double eval(sgpp::base::DataVector const& v) { return v[0]; }

}  // namespace debugfct
// ----------------------------------------------------------------------------------

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
  double tolerance = 2e-8;
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

// This is only used for Debug purposes.
// The variance data is read from a file to make it easily replacable
BOOST_AUTO_TEST_CASE(testVarianceDEBUG) {
  std::cout << "Testing B spline variance calculation  subgridwise on single levels from file"
            << std::endl;
  std::ifstream varianceFile("combigrid/tests/BSplineVarianceData.dat");
  std::string line;
  size_t level1;
  size_t level2;
  double variance;
  std::vector<double> varianceError;
  std::vector<sgpp::combigrid::MultiIndex> levels;
  while (varianceFile >> level1 >> level2 >> variance) {
    sgpp::combigrid::MultiIndex level{level1, level2};
    double bSplineVariance = BSplineVariance(level);
    varianceError.push_back(std::fabs(bSplineVariance - variance));
    levels.push_back(level);
  }
  for (size_t i = 0; i < levels.size(); i++) {
    std::cout << levels[i][0] << " " << levels[i][1] << "|  error:  " << varianceError[i]
              << std::endl;
  }
  varianceFile.close();
}

BOOST_AUTO_TEST_CASE(testVarianceOnDiagonal) {
  std::cout
      << "Testing B spline variance calculation on levels of the diagonal of the subgrid scheme"
      << std::endl;
  VarianceDiagonalTestData varianceDiagonalTestData;
  double varianceError;

  for (size_t i = 0; i < varianceDiagonalTestData.levels.size(); i++) {
    MultiIndex level = varianceDiagonalTestData.levels[i];
    double bSplineVariance = BSplineVariance(level);
    varianceError = std::fabs(bSplineVariance - varianceDiagonalTestData.variance);
    std::cout << "level: " << level[0] << " " << level[1] << "|  error:  " << varianceError
              << std::endl;
    //    BOOST_CHECK_SMALL(varianceError, 1e-6);
  }
}

/*
 * This test creates a Bspline interpolant and calculates mean and variance with Monte Carlo methods
 * for two test functions for which mean and variance are known
 *
 * Comment: The results obtained by Monte Carlo are extremely unprecise
 */
BOOST_AUTO_TEST_CASE(testMCVariance) {
  std::cout << "Testing Variance of B spline interpolants using Monte Carlo methods" << std::endl;
  size_t numDimensions = debugfct::numDims;
  size_t degree = 3;
  sgpp::combigrid::MultiFunction func(debugfct::eval);

  std::vector<std::shared_ptr<sgpp::combigrid::AbstractPointHierarchy>> pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  std::vector<
      std::shared_ptr<sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatArrayVector>>>
      evaluators(numDimensions,
                 sgpp::combigrid::CombiEvaluators::multiBSplineInterpolation(degree));
  //  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
  //      new sgpp::combigrid::AveragingLevelManager());
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::WeightedRatioLevelManager());
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);
  bool exploitNesting = false;

  auto multiOperation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, evaluators, levelManager, gf, exploitNesting);
  double diff = 0.0;
  // generator generates num_points random points in [0,1]^numDimensions
  size_t num_MCpoints = 10000;
  sgpp::quadrature::NaiveSampleGenerator generator(numDimensions);
  sgpp::base::DataVector p(numDimensions, 0);
  std::vector<sgpp::base::DataVector> params;

  for (size_t i = 0; i < num_MCpoints; i++) {
    generator.getSample(p);
    params.push_back(p);
  }

  for (size_t level = 0; level < 11; level++) {
    std::cout << "level " << level << ": ";
    double MCL2_err = 0.0;
    double MCMean = 0.0;
    double MCVar = 0.0;
    auto result = multiOperation->evaluate(level, params);
    for (size_t i = 0; i < params.size(); ++i) {
      diff = func(params[i]) - result[i];
      MCL2_err += diff * diff;
      MCMean += result[i];
      MCVar += std::pow(result[i] - debugfct::mean, 2);
    }
    MCMean /= static_cast<double>(num_MCpoints);
    MCVar /= static_cast<double>(num_MCpoints);

    MCL2_err = sqrt(MCL2_err / static_cast<double>(num_MCpoints));
    double MCMean_err = fabs(MCMean - debugfct::mean);
    double MCVar_err = fabs(MCVar - debugfct::variance);

    //    std::cout << MCL2_err << " " << MCMean_err << " " << MCVar_err << std::endl;
    std::cout << MCL2_err << " " << MCMean_err << " " << MCVar_err << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
