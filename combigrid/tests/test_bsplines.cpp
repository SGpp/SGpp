// Copyright(C)2008 - today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridGridBasedEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridCallbackEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/quadrature/sampling/NaiveSampleGenerator.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

double BSplineVariance(sgpp::combigrid::MultiIndex level) {
  sgpp::combigrid::Atan atanModel;
  size_t numDimensions = atanModel.numDims;
  size_t degree = 3;
  sgpp::combigrid::MultiFunction func(atanModel.eval);
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

double polynomialVarianceQuadrature(sgpp::combigrid::MultiIndex& level) {
  sgpp::combigrid::Atan atanModel;
  size_t numDimensions = atanModel.numDims;
  sgpp::combigrid::MultiFunction func(atanModel.eval);
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expClenshawCurtis());

  sgpp::combigrid::EvaluatorConfiguration evalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_PolynomialScalarProduct);

  sgpp::combigrid::CombiEvaluators::MultiCollection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(evalConfig));

  bool exploitNesting = true;
  auto summationStrategyType = sgpp::combigrid::FullGridSummationStrategyType::VARIANCE;

  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> storage =
      std::make_shared<sgpp::combigrid::CombigridTreeStorage>(pointHierarchies, exploitNesting,
                                                              func);

  auto fullGridEval = std::make_shared<
      sgpp::combigrid::FullGridCallbackEvaluator<sgpp::combigrid::FloatArrayVector>>(
      storage, evaluators, pointHierarchies, summationStrategyType);

  auto result = fullGridEval->eval(level);
  double res = result[0].value();
  return res;
}

double polynomialVariancePCE(sgpp::combigrid::MultiIndex& level) {
  sgpp::combigrid::Atan atanModel;
  size_t numDimensions = atanModel.numDims;
  sgpp::combigrid::MultiFunction func(atanModel.eval);
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expClenshawCurtis());

  sgpp::combigrid::EvaluatorConfiguration evalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Tensor_PolynomialInterpolation);
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  evalConfig.functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::CombiEvaluators::TensorCollection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiTensorEvaluator(evalConfig));

  auto summationStrategyType = sgpp::combigrid::FullGridSummationStrategyType::LINEAR;

  auto storage =
      std::make_shared<sgpp::combigrid::CombigridTreeStorage>(pointHierarchies, true, func);

  auto fullGridEval = std::make_shared<
      sgpp::combigrid::FullGridCallbackEvaluator<sgpp::combigrid::FloatTensorVector>>(
      storage, evaluators, pointHierarchies, summationStrategyType);

  return fullGridEval->eval(level).norm();
}

BOOST_AUTO_TEST_SUITE(testPolynomialVariance)

BOOST_AUTO_TEST_CASE(testVarianceOfPolynomialsOnDiagonal) {
  sgpp::combigrid::Atan atanModel;
  std::vector<double> tolerancesQuad{1e1, 1e0, 1e0, 1e0, 1e-1, 1e-2, 1e-3, 1e-5};
  std::vector<double> tolerancesPCE{1e1, 1e0, 1e0, 1e0, 1e-1, 1e-2, 1e-3, 1e-5};

  sgpp::combigrid::MultiIndex level(atanModel.numDims);
  for (size_t i = 0; i < 7; i++) {
    sgpp::combigrid::MultiIndex level(atanModel.numDims, i);
    double polyVariance = polynomialVarianceQuadrature(level);
    double varianceError = std::fabs(polyVariance - atanModel.variance);
    //    std::cout << "level: |" << level[0] << " " << level[1] << "| error:  " << varianceError
    //              << std::endl;
    BOOST_CHECK_SMALL(varianceError, tolerancesQuad[i]);

#ifdef USE_DAKOTA
    // -----------------------------------------------------------------------
    // use the PCE approach
    polyVariance = polynomialVariancePCE(level);
    varianceError = std::fabs(polyVariance - atanModel.variance);
    //    std::cout << "level: |" << level[0] << " " << level[1] << "| error:  " << varianceError
    //              << std::endl;
    BOOST_CHECK_SMALL(varianceError, tolerancesPCE[i]);
  }
#endif
}

BOOST_AUTO_TEST_CASE(testVarianceOfPolynomialsOnLevel) {
  // This data was created by SGpp/combigrid/tests/createVarianceDataPolynomial.py
  // It represents the variance calculated for each level with levelsum <= 7 for the function
  // arctan(50.0 * (x[0] - .35)) + pi / 2.0 + 4.0 * x[1] ** 3 + exp(x[0] * x[1] - 1.0)
  struct AtanModelVarianceTestDataPolynomials {
    std::vector<sgpp::combigrid::MultiIndex> levels{
        sgpp::combigrid::MultiIndex{1, 3}, sgpp::combigrid::MultiIndex{3, 0},
        sgpp::combigrid::MultiIndex{0, 7}, sgpp::combigrid::MultiIndex{1, 6},
        sgpp::combigrid::MultiIndex{5, 1}, sgpp::combigrid::MultiIndex{2, 5},
        sgpp::combigrid::MultiIndex{0, 3}, sgpp::combigrid::MultiIndex{4, 0},
        sgpp::combigrid::MultiIndex{1, 2}, sgpp::combigrid::MultiIndex{3, 3},
        sgpp::combigrid::MultiIndex{2, 0}, sgpp::combigrid::MultiIndex{1, 5},
        sgpp::combigrid::MultiIndex{5, 0}, sgpp::combigrid::MultiIndex{2, 2},
        sgpp::combigrid::MultiIndex{4, 1}, sgpp::combigrid::MultiIndex{1, 1},
        sgpp::combigrid::MultiIndex{3, 2}, sgpp::combigrid::MultiIndex{0, 0},
        sgpp::combigrid::MultiIndex{0, 4}, sgpp::combigrid::MultiIndex{6, 0},
        sgpp::combigrid::MultiIndex{1, 4}, sgpp::combigrid::MultiIndex{2, 3},
        sgpp::combigrid::MultiIndex{2, 1}, sgpp::combigrid::MultiIndex{4, 2},
        sgpp::combigrid::MultiIndex{1, 0}, sgpp::combigrid::MultiIndex{0, 1},
        sgpp::combigrid::MultiIndex{7, 0}, sgpp::combigrid::MultiIndex{5, 2},
        sgpp::combigrid::MultiIndex{6, 1}, sgpp::combigrid::MultiIndex{3, 1},
        sgpp::combigrid::MultiIndex{0, 2}, sgpp::combigrid::MultiIndex{0, 6},
        sgpp::combigrid::MultiIndex{4, 3}, sgpp::combigrid::MultiIndex{0, 5},
        sgpp::combigrid::MultiIndex{3, 4}, sgpp::combigrid::MultiIndex{2, 4}};
    std::vector<double> variances{
        2.549880259369381, 2.071801193244269,  1.437020666045363, 2.549880259369484,
        3.711520362783947, 3.343415383888980,  1.437020666045363, 2.020867648161442,
        2.549879363026577, 3.545669384171749,  1.871362920523577, 2.549880259369488,
        1.971865046276090, 3.343414644093775,  3.760627367537953, 2.816790782428100,
        3.545668634996716, -0.000000000000004, 1.437020666045356, 1.980062220760152,
        2.549880259369484, 3.343415383888923,  3.610252425086911, 3.493781280710625,
        1.080108355131403, 1.701156850402693,  1.980270791279654, 3.444674880419520,
        3.719739391135988, 3.812521444326839,  1.437020562734489, 1.437020666045367,
        3.493782026820650, 1.437020666045360,  3.545669384171807, 3.343415383888976};
  };

  AtanModelVarianceTestDataPolynomials varianceTestData;
  for (size_t i = 0; i < varianceTestData.levels.size(); i++) {
    sgpp::combigrid::MultiIndex level = varianceTestData.levels[i];

    // -----------------------------------------------------------------------
    // use the quadrature approach
    double polyVariance = polynomialVarianceQuadrature(level);
    double varianceError = std::fabs(polyVariance - varianceTestData.variances[i]);
    //    std::cout << "level: |" << level[0] << " " << level[1] << "|  error:  " << varianceError
    //              << std::endl;
    BOOST_CHECK_SMALL(varianceError, 5e-13);

#ifdef USE_DAKOTA
    // -----------------------------------------------------------------------
    // use the PCE approach
    polyVariance = polynomialVariancePCE(level);
    varianceError = std::fabs(polyVariance - varianceTestData.variances[i]);
    //    std::cout << "level: |" << level[0] << " " << level[1] << "|  value: " << polyVariance
    //              << " (err=" << varianceError << ")" << std::endl;
    BOOST_CHECK_SMALL(varianceError, 1e-14);
  }
#endif
}

BOOST_AUTO_TEST_SUITE_END()

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

  // This data was created by SGpp/combigrid/tests/createVarianceDataBsplines.py
  // It represents the variance calculated for each level with levelsum <= 7 for the function
  // arctan(50.0 * (x[0] - .35)) + pi / 2.0 + 4.0 * x[1] ** 3 + exp(x[0] * x[1] - 1.0)
  struct AtanModelVarianceTestDataBsplines {
    std::vector<sgpp::combigrid::MultiIndex> levels{
        sgpp::combigrid::MultiIndex{0, 0}, sgpp::combigrid::MultiIndex{0, 1},
        sgpp::combigrid::MultiIndex{0, 2}, sgpp::combigrid::MultiIndex{0, 3},
        sgpp::combigrid::MultiIndex{0, 4}, sgpp::combigrid::MultiIndex{0, 5},
        sgpp::combigrid::MultiIndex{0, 6}, sgpp::combigrid::MultiIndex{0, 7},
        sgpp::combigrid::MultiIndex{1, 0}, sgpp::combigrid::MultiIndex{1, 1},
        sgpp::combigrid::MultiIndex{1, 2}, sgpp::combigrid::MultiIndex{1, 3},
        sgpp::combigrid::MultiIndex{1, 4}, sgpp::combigrid::MultiIndex{1, 5},
        sgpp::combigrid::MultiIndex{1, 6}, sgpp::combigrid::MultiIndex{2, 0},
        sgpp::combigrid::MultiIndex{2, 1}, sgpp::combigrid::MultiIndex{2, 2},
        sgpp::combigrid::MultiIndex{2, 3}, sgpp::combigrid::MultiIndex{2, 4},
        sgpp::combigrid::MultiIndex{2, 5}, sgpp::combigrid::MultiIndex{3, 0},
        sgpp::combigrid::MultiIndex{3, 1}, sgpp::combigrid::MultiIndex{3, 2},
        sgpp::combigrid::MultiIndex{3, 3}, sgpp::combigrid::MultiIndex{3, 4},
        sgpp::combigrid::MultiIndex{4, 0}, sgpp::combigrid::MultiIndex{4, 1},
        sgpp::combigrid::MultiIndex{4, 2}, sgpp::combigrid::MultiIndex{4, 3},
        sgpp::combigrid::MultiIndex{5, 0}, sgpp::combigrid::MultiIndex{5, 1},
        sgpp::combigrid::MultiIndex{5, 2}, sgpp::combigrid::MultiIndex{6, 0},
        sgpp::combigrid::MultiIndex{6, 1}, sgpp::combigrid::MultiIndex{7, 0}};
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

  AtanModelVarianceTestDataBsplines varianceTestData;
  double tolerance = 2e-8;
  for (size_t i = 0; i < varianceTestData.levels.size(); i++) {
    sgpp::combigrid::MultiIndex level = varianceTestData.levels[i];
    double bSplineVariance = BSplineVariance(level);
    double varianceError = std::fabs(bSplineVariance - varianceTestData.variances[i]);
    std::cout << "level: " << level[0] << " " << level[1] << "|  error:  " << varianceError
              << std::endl;
    //    BOOST_CHECK_SMALL(varianceError, tolerance);
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
  sgpp::combigrid::Atan atanModel;
  for (size_t i = 0; i < 7; i++) {
    sgpp::combigrid::MultiIndex level(atanModel.numDims, i);
    double bSplineVariance = BSplineVariance(level);
    double varianceError = std::fabs(bSplineVariance - atanModel.variance);
    std::cout << "level: " << level[0] << " " << level[1] << "|  error:  " << varianceError
              << std::endl;
    //    BOOST_CHECK_SMALL(varianceError, 1e-6);
  }
}

/*
 * This test creates a Bspline interpolant and calculates mean and variance with Monte Carlo methods
 * for two test functions for which mean and variance are known
 *
 * Comment: The results obtained by Monte Carlo are extremely unprecise even for simple functions as
 * f(x)=x
 */
BOOST_AUTO_TEST_CASE(testMCVariance) {
  std::cout << "Testing Variance of B spline interpolants using Monte Carlo methods" << std::endl;
  sgpp::combigrid::Debugfct debugModel;
  size_t numDimensions = debugModel.numDims;
  size_t degree = 3;
  sgpp::combigrid::MultiFunction func(debugModel.eval);

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
      MCVar += std::pow(result[i] - debugModel.mean, 2);
    }
    MCMean /= static_cast<double>(num_MCpoints);
    MCVar /= static_cast<double>(num_MCpoints);

    MCL2_err = sqrt(MCL2_err / static_cast<double>(num_MCpoints));
    double MCMean_err = fabs(MCMean - debugModel.mean);
    double MCVar_err = fabs(MCVar - debugModel.variance);

    //    std::cout << MCL2_err << " " << MCMean_err << " " << MCVar_err << std::endl;
    std::cout << MCL2_err << " " << MCMean_err << " " << MCVar_err << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
