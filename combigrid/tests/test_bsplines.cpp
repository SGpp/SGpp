// Copyright(C)2008 - today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridCallbackEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridGridBasedEvaluator.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotNakBsplineBoundaryCombigrid.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/quadrature/sampling/NaiveSampleGenerator.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

double BSplineVariance(sgpp::combigrid::MultiIndex level, size_t degree) {
  sgpp::combigrid::AtanUniform atanModel;
  size_t numDimensions = atanModel.numDims;
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
  sgpp::combigrid::AtanUniform atanModel;
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
  sgpp::combigrid::AtanUniform atanModel;
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

  // compute the variance
  auto tensor = fullGridEval->eval(level);
  double mean = tensor.get(sgpp::combigrid::MultiIndex(numDimensions, 0)).value();
  double variance = std::pow(tensor.norm(), 2) - std::pow(mean, 2);

  return variance;
}

BOOST_AUTO_TEST_SUITE(testPolynomialVariance)

BOOST_AUTO_TEST_CASE(testVarianceOfPolynomialsOnDiagonal) {
  sgpp::combigrid::AtanUniform atanModel;
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
#endif
  }
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
    BOOST_CHECK_SMALL(varianceError, 5e-14);
#endif
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(testBsplines)

double x1(sgpp::base::DataVector const& v) { return v[0] + v[1]; }
double x3(sgpp::base::DataVector const& v) { return std::pow(v[0], 3) + std::pow(v[1], 3); }
double x5(sgpp::base::DataVector const& v) { return std::pow(v[0], 5) + std::pow(v[1], 5); }

double L2BsplineInterpolationError(size_t numDimensions, size_t degree,
                                   sgpp::combigrid::MultiFunction func, size_t level) {
  auto operation =
      sgpp::combigrid::CombigridMultiOperation::createExpUniformBoundaryBsplineInterpolation(
          numDimensions, func, degree);

  double L2Err = 0.0;
  size_t numMCpoints = 1000;
  sgpp::quadrature::NaiveSampleGenerator generator(numDimensions);
  sgpp::base::DataMatrix params(numDimensions, numMCpoints);
  sgpp::base::DataVector p(numDimensions);
  sgpp::base::DataVector Feval(numMCpoints, 0.0);
  for (size_t i = 0; i < numMCpoints; i++) {
    generator.getSample(p);
    params.setColumn(i, p);
    Feval.set(i, func(p));
  }
  operation->setParameters(params);
  operation->getLevelManager()->addRegularLevels(level);
  sgpp::base::DataVector eval = operation->getResult();
  eval.sub(Feval);
  for (size_t i = 0; i < eval.size(); i++) {
    L2Err += fabs(eval[i] * eval[i]);
  }
  L2Err = sqrt(L2Err / static_cast<double>(numMCpoints));
  return L2Err;
}

BOOST_AUTO_TEST_CASE(testCorrespondingDegreeInterpolation) {
  std::cout
      << "Testing interpolation of x^d+y^d for B splines of degree d on level d, d in {1,3,5}."
      << std::endl;
  size_t numDimensions = 2;
  size_t degree = 1;
  sgpp::combigrid::MultiFunction func1(x1);
  size_t level = 1;
  double L2error1 = L2BsplineInterpolationError(numDimensions, degree, func1, level);
  degree = 3;
  sgpp::combigrid::MultiFunction func3(x3);
  level = 3;
  double L2error3 = L2BsplineInterpolationError(numDimensions, degree, func3, level);
  degree = 5;
  sgpp::combigrid::MultiFunction func5(x5);
  level = 5;
  double L2error5 = L2BsplineInterpolationError(numDimensions, degree, func5, level);
  double tolerance = 3e-14;

  //  std::cout << "d = 1: " << L2error1 << std::endl;
  //  std::cout << "d = 3: " << L2error3 << std::endl;
  //  std::cout << "d = 5: " << L2error5 << "\n" << std::endl;

  BOOST_CHECK_SMALL(L2error1, tolerance);
  BOOST_CHECK_SMALL(L2error3, tolerance);
  BOOST_CHECK_SMALL(L2error5, tolerance);
}

double BsplineQuadratureError(size_t numDimensions, size_t degree,
                              sgpp::combigrid::MultiFunction func, size_t level) {
  auto operation = sgpp::combigrid::CombigridOperation::createExpUniformBoundaryBsplineQuadrature(
      numDimensions, func, degree);
  double integral = operation->evaluate(level);
  // correct solution is 2/(degree+1)
  return fabs(integral - 2.0 / (static_cast<double>(degree) + 1));
}
BOOST_AUTO_TEST_CASE(testCorrespondingDegreeQuadrature) {
  std::cout << "Testing integration of x^d+y^d for B splines of degree d on level d, d in {1,3,5}."
            << std::endl;
  size_t numDimensions = 2;
  size_t degree = 1;
  sgpp::combigrid::MultiFunction func1(x1);
  size_t level = 1;
  double Quaderror1 = BsplineQuadratureError(numDimensions, degree, func1, level);
  degree = 3;
  sgpp::combigrid::MultiFunction func3(x3);
  level = 3;
  double Quaderror3 = BsplineQuadratureError(numDimensions, degree, func3, level);
  degree = 5;
  sgpp::combigrid::MultiFunction func5(x5);
  level = 5;
  double Quaderror5 = BsplineQuadratureError(numDimensions, degree, func5, level);
  double tolerance = 1e-14;

  //  std::cout << "d = 1: " << Quaderror1 << std::endl;
  //  std::cout << "d = 3: " << Quaderror3 << std::endl;
  //  std::cout << "d = 5: " << Quaderror5 << "\n" << std::endl;

  BOOST_CHECK_SMALL(Quaderror1, tolerance);
  BOOST_CHECK_SMALL(Quaderror3, tolerance);
  BOOST_CHECK_SMALL(Quaderror5, tolerance);
}

// Does this test belong here? (OperationMatrixLTwoDotNakBsplineBoundaryCombigrid is defined in pde
// module but only used for combigrids)

double BsplineQuadratureSquare(size_t numDimensions, size_t degree,
                               sgpp::combigrid::MultiFunction func, size_t level) {
  std::shared_ptr<sgpp::base::Grid> grid;
  grid.reset(sgpp::base::Grid::createNakBsplineBoundaryCombigridGrid(numDimensions, degree));
  grid->getGenerator().regular(level);
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::optimization::sle_solver::Auto sleSolver;
  sgpp::base::DataVector alpha(grid->getSize());
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  sgpp::base::DataVector f_values(gridStorage.getSize(), 0.0);
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
      p[j] = gp.getStandardCoordinate(j);
    }
    f_values[i] = func(p);
  }

  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
    std::cout << "Solving failed!" << std::endl;
  }
  sgpp::base::Grid* gridptr = grid.get();
  sgpp::pde::OperationMatrixLTwoDotNakBsplineBoundaryCombigrid massMatrix(gridptr);
  sgpp::base::DataVector product(alpha.size(), 0);
  massMatrix.mult(alpha, product);
  double integralSquare = product.dotProduct(alpha);
  return integralSquare;
}
BOOST_AUTO_TEST_CASE(testCorrespondingDegreeScalarProducts) {
  std::cout
      << "Testing integration of (x^d+y^d)^2 for B splines of degree d on level d, d in {1,3,5}.\n"
      << "This verifies the scalar product routine." << std::endl;
  size_t numDimensions = 2;
  size_t degree = 1;
  sgpp::combigrid::MultiFunction func1(x1);
  size_t level = 1;
  double QuadSquare1 = BsplineQuadratureSquare(numDimensions, degree, func1, level);
  double QuadSquareError1 = fabs(QuadSquare1 - 7.0 / 6.0);
  degree = 3;
  sgpp::combigrid::MultiFunction func3(x3);
  level = 3;
  double QuadSquare3 = BsplineQuadratureSquare(numDimensions, degree, func3, level);
  double QuadSquareError3 = fabs(QuadSquare3 - 23.0 / 56.0);
  degree = 5;
  sgpp::combigrid::MultiFunction func5(x5);
  level = 5;
  double QuadSquare5 = BsplineQuadratureSquare(numDimensions, degree, func5, level);
  double QuadSquareError5 = fabs(QuadSquare5 - 47.0 / 198.0);
  double tolerance = 1e-14;

  //  std::cout << "d = 1: " << QuadSquareError1 << std::endl;
  //  std::cout << "d = 3: " << QuadSquareError3 << std::endl;
  //  std::cout << "d = 5: " << QuadSquareError5 << "\n" << std::endl;

  BOOST_CHECK_SMALL(QuadSquareError1, tolerance);
  BOOST_CHECK_SMALL(QuadSquareError3, tolerance);
  BOOST_CHECK_SMALL(QuadSquareError5, tolerance);
}

/*
 * This test calculates the variance for the test function atanModel (see above) and compares it
 * to precalculated data.
 *
 */

BOOST_AUTO_TEST_CASE(testVarianceOnLeveldeg3) {
  std::cout << "Testing B spline variance calculation  subgridwise on single levels for B splines "
               "of degree 3."
            << std::endl;

  // This data was created by SGpp/combigrid/tests/createVarianceDataBsplines.py
  // It represents the variance calculated for each level with levelsum <= 7 for the function
  // arctan(50.0 * (x[0] - .35)) + pi / 2.0 + 4.0 * x[1] ** 3 + exp(x[0] * x[1] - 1.0)
  struct AtanModelVarianceTestDataBsplines {
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
        2.549880760493220, 1.858681146494114,  1.437020666045346, 2.549880259378867,
        3.719658912837517, 3.900926591786073,  1.437020728026305, 1.985491657343513,
        2.549893635809379, 3.330801307994530,  2.427023144889512, 2.549880259715540,
        1.979975007603208, 3.900939092268690,  3.725129287185531, 2.816790782428100,
        3.330812831246485, -0.000000000000004, 1.437020668119384, 1.980262657807472,
        2.549880271385291, 3.900927051314765,  4.167777616534886, 3.458296945247319,
        1.080108355131403, 1.701156850402690,  1.980269030387413, 3.452826245565596,
        3.719940998350237, 3.597641299962191,  1.437022054119371, 1.437020666048241,
        3.458285162235061, 1.437020666114215,  3.330800874338360, 3.900926602175467};
  };

  AtanModelVarianceTestDataBsplines varianceTestData;
  double tolerance = 2e-8;
  for (size_t i = 0; i < varianceTestData.levels.size(); i++) {
    sgpp::combigrid::MultiIndex level = varianceTestData.levels[i];
    size_t degree = 3;
    double bSplineVariance = BSplineVariance(level, degree);
    double varianceError = std::fabs(bSplineVariance - varianceTestData.variances[i]);
    //    std::cout << "level: " << level[0] << " " << level[1] << "|  error:  " << varianceError
    //              << std::endl;
    BOOST_CHECK_SMALL(varianceError, tolerance);
  }
}

BOOST_AUTO_TEST_CASE(testVarianceOnLeveldeg5) {
  std::cout << "Testing B spline variance calculation  subgridwise on single levels for B splines "
               "of degree 5."
            << std::endl;

  // This data was created by SGpp/combigrid/tests/createVarianceDataBsplines.py
  // It represents the variance calculated for each level with levelsum <= 7 for the function
  // arctan(50.0 * (x[0] - .35)) + pi / 2.0 + 4.0 * x[1] ** 3 + exp(x[0] * x[1] - 1.0)
  struct AtanModelVarianceTestDataBsplines {
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
        2.549880263014657, 1.758709241121634,  1.437020666045367, 2.549880259369480,
        3.719898759162344, 4.019259365726171,  1.437020666146449, 2.001767995410288,
        2.549883454168121, 3.230298971637547,  2.545726212447136, 2.549880259369768,
        1.980214864795903, 4.019262025441856,  3.741453545290554, 2.816790782428100,
        3.230301616151564, -0.000000000000004, 1.437020666046340, 1.980283267070979,
        2.549880259403881, 4.019259368387642,  4.286103665868156, 3.474611382611405,
        1.080108355131403, 1.701156850402690,  1.980270151275377, 3.453056488503831,
        3.719961616942067, 3.497136996312333,  1.437021035230856, 1.437020666045367,
        3.474608731181096, 1.437020666045370,  3.230298969095825, 4.019259365750923};
  };

  AtanModelVarianceTestDataBsplines varianceTestData;
  double tolerance = 1e-8;
  for (size_t i = 0; i < varianceTestData.levels.size(); i++) {
    sgpp::combigrid::MultiIndex level = varianceTestData.levels[i];
    size_t degree = 5;
    double bSplineVariance = BSplineVariance(level, degree);
    double varianceError = std::fabs(bSplineVariance - varianceTestData.variances[i]);
    //    std::cout << "level: " << level[0] << " " << level[1] << "|  error:  " << varianceError
    //              << std::endl;
    BOOST_CHECK_SMALL(varianceError, tolerance);
  }
}

BOOST_AUTO_TEST_CASE(testVarianceOnDiagonaldeg3) {
  std::cout << "Testing B spline variance calculation on levels of the diagonal of the subgrid "
               "scheme for degree 3."
            << std::endl;
  sgpp::combigrid::AtanUniform atanModel;
  std::vector<double> tolerance = {4, 1, 0.5, 0.2, 0.006, 0.0003, 8e-06};
  size_t degree = 3;
  for (size_t i = 0; i < 7; i++) {
    sgpp::combigrid::MultiIndex level(atanModel.numDims, i);
    double bSplineVariance = BSplineVariance(level, degree);
    double varianceError = std::fabs(bSplineVariance - atanModel.variance);
    //    std::cout << "level: " << level[0] << " " << level[1] << "|  error:  " << varianceError
    //              << std::endl;
    BOOST_CHECK_SMALL(varianceError, tolerance[i]);
  }
}

BOOST_AUTO_TEST_CASE(testVarianceOnDiagonaldeg5) {
  std::cout << "Testing B spline variance calculation on levels of the diagonal of the subgrid "
               "scheme for degree 5."
            << std::endl;
  sgpp::combigrid::AtanUniform atanModel;
  std::vector<double> tolerance = {4, 1, 0.6, 0.25, 0.022, 5e-05, 2e-05};
  size_t degree = 5;
  for (size_t i = 0; i < 7; i++) {
    sgpp::combigrid::MultiIndex level(atanModel.numDims, i);
    double bSplineVariance = BSplineVariance(level, degree);
    double varianceError = std::fabs(bSplineVariance - atanModel.variance);
    //    std::cout << "level: " << level[0] << " " << level[1] << "|  error:  " << varianceError
    //              << std::endl;
    BOOST_CHECK_SMALL(varianceError, tolerance[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
