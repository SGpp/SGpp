// Copyright(C)2008 - today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/functions/ProbabilityDensityFunction1D.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/RegularLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridCallbackEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridGridBasedEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/sparsegrid/LTwoScalarProductHashMapNakBsplineBoundaryCombigrid.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineScalarProductEvaluator.hpp>
#include <sgpp/combigrid/pce/BsplineStochasticCollocation.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>
#include "../../base/src/sgpp/base/tools/sle/solver/Auto.hpp"
#include "../../base/src/sgpp/base/tools/sle/system/HierarchisationSLE.hpp"

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

// auxiliary function. Creates a regular level structure and calculated the Bspline interpolation
// coefficients for this structure
void createRegularLevelStructure(
    size_t numLevels, size_t degree,
    sgpp::combigrid::CombiHierarchies::Collection const& pointHierarchies,
    sgpp::combigrid::GridFunction gf, bool exploitNesting,
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>>& newLevelStructure,
    std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage>& coefficientStorage,
    size_t numDimensions) {
  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, degree);
  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::RegularLevelManager());
  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;

  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, Evaluators, levelManager, gf, exploitNesting,
      auxiliarySummationStrategyType);
  Operation->getLevelManager()->addRegularLevels(numLevels);
  coefficientStorage = Operation->getStorage();
  newLevelStructure = Operation->getLevelManager()->getLevelStructure();
}

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

BOOST_AUTO_TEST_CASE(testBSCInterpolation) {
  std::cout << "Testing Bspline Stochastic Collocation's interpolation eval." << std::endl;

  size_t numDimensions = 2;
  sgpp::combigrid::MultiFunction func(x1);
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> storage;
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager =
      std::make_shared<sgpp::combigrid::AveragingLevelManager>();
  sgpp::combigrid::CombigridSurrogateModelConfiguration bsc_config;
  bsc_config.type = sgpp::combigrid::CombigridSurrogateModelsType::BSPLINE_STOCHASTIC_COLLOCATION;
  bsc_config.pointHierarchies = pointHierarchies;
  bsc_config.levelManager = levelManager;

  bsc_config.degree = 5;
  bsc_config.numDimensions = numDimensions;
  bsc_config.coefficientStorage = storage;

  sgpp::combigrid::BsplineStochasticCollocation bsc(bsc_config);

  // create level Structure and interpolate
  std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> newLevelStructure;
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> newCoefficientStorage;
  size_t numLevels = 5;
  sgpp::combigrid::GridFunction gf =
      BSplineCoefficientGridFunction(func, pointHierarchies, bsc_config.degree);
  createRegularLevelStructure(numLevels, bsc_config.degree, pointHierarchies, gf, false,
                              newLevelStructure, newCoefficientStorage, numDimensions);

  //  update config
  bsc_config.levelStructure = newLevelStructure;
  bsc_config.coefficientStorage = newCoefficientStorage;

  bsc.updateConfig(bsc_config);

  // 0.1  0.2  0.3
  // 0.4  0.5  0.6
  sgpp::base::DataMatrix pointsT(2, 3);
  pointsT.set(0, 0, 0.1);
  pointsT.set(0, 1, 0.2);
  pointsT.set(0, 2, 0.3);
  pointsT.set(1, 0, 0.4);
  pointsT.set(1, 1, 0.5);
  pointsT.set(1, 2, 0.6);
  sgpp::base::DataVector res(3);
  bsc.eval(pointsT, res);

  //  std::cout << res[0] << " " << res[1] << " " << res[2] << std::endl;
  double tolerance = 1e-13;
  BOOST_CHECK_SMALL(res[0] - 0.5, tolerance);
  BOOST_CHECK_SMALL(res[1] - 0.7, tolerance);
  BOOST_CHECK_SMALL(res[2] - 0.9, tolerance);

  sgpp::base::DataVector point;
  point.append(0.337);
  point.append(0.9919);
  double res1 = bsc.eval(point);
  //  std::cout << res1 << std::endl;
  BOOST_CHECK_SMALL(res1 - 1.3289, tolerance);
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

  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
    std::cout << "Solving failed!" << std::endl;
  }
  sgpp::base::Grid* gridptr = grid.get();
  sgpp::combigrid::LTwoScalarProductHashMapNakBsplineBoundaryCombigrid massMatrix(gridptr);
  sgpp::base::DataVector product(alpha.size(), 0);
  massMatrix.mult(alpha, product);
  double integralSquare = product.dotProduct(alpha);
  return integralSquare;
}
BOOST_AUTO_TEST_CASE(testCorrespondingDegreeOperationMatrixScalarProducts) {
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

double omscFunc(sgpp::base::DataVector const& v) { return std::pow(v[0], 5) + std::pow(v[1], 5); }
double omscWeight(double x) { return sin(2 * x); }  // <= sin(x) on [0,2]^2

BOOST_AUTO_TEST_CASE(testBsplineLTwoScalarProductsWithWeightsAndBounds) {
  std::cout << "calculating mean and variance of x^5 + y^5 on [0,1]^2 using Bspline scalar "
               "products and the weight function sin(x) on [0,2]^2."
            << std::endl;
  size_t numDims = 2;
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDims, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> storage;
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::RegularLevelManager());
  sgpp::combigrid::SingleFunction weightfunction(omscWeight);
  sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection(0);
  sgpp::base::DataVector bounds(0);
  for (size_t d = 0; d < numDims; d++) {
    bounds.push_back(0);
    bounds.push_back(2);
    weightFunctionsCollection.push_back(weightfunction);
  }

  sgpp::combigrid::CombigridSurrogateModelConfiguration config;
  config.type = sgpp::combigrid::CombigridSurrogateModelsType::BSPLINE_STOCHASTIC_COLLOCATION;
  config.pointHierarchies = pointHierarchies;
  config.storage = storage;
  config.levelManager = levelManager;
  config.degree = 3;
  config.coefficientStorage = storage;
  config.weightFunctions = weightFunctionsCollection;

  config.bounds = bounds;
  sgpp::combigrid::BsplineStochasticCollocation BSC(config);

  std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> levelStructure;
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> coefficientStorage;
  sgpp::combigrid::MultiFunction func(omscFunc);
  sgpp::combigrid::GridFunction gf =
      BSplineCoefficientGridFunction(func, pointHierarchies, config.degree);

  size_t numLevels = 5;

  createRegularLevelStructure(numLevels, config.degree, pointHierarchies, gf, false, levelStructure,
                              coefficientStorage, numDims);
  //  update config
  config.levelStructure = levelStructure;
  config.coefficientStorage = coefficientStorage;
  BSC.updateConfig(config);

  // x^5+y^5 on [0,1]^2 with weight function sin(x) on [0,2]^2
  double realEv = 0.906028496608237;
  double realVar = 0.700571273115382;

  double var = BSC.variance();
  double ev = BSC.mean();
  //  std::cout << fabs(ev - realEv) << " " << fabs(var - realVar) << std::endl;
  BOOST_CHECK_SMALL(fabs(ev - realEv), 1e-6);
  BOOST_CHECK_SMALL(fabs(var - realVar), 1e-6);
}

/*
 * This test calculates the variance for the test function atanModel  and compares it
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

double x32D(sgpp::base::DataVector const& v) { return std::pow(v[0], 3) + std::pow(v[1], 3); }
double wcos(double v) { return cos(v); }

BOOST_AUTO_TEST_CASE(testQuadratureWithWeightFunction) {
  std::cout << "Integrate objective function x^3+y^3 and weight function cos x with B splines of "
               "degree 3 on level 5 in 2D "
            << std::endl;
  size_t numDimensions = 2;
  size_t degree = 3;
  sgpp::combigrid::MultiFunction func(x32D);
  sgpp::combigrid::SingleFunction weightfunction(wcos);
  sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection(0);
  for (size_t d = 0; d < numDimensions; d++) {
    weightFunctionsCollection.push_back(weightfunction);
  }
  size_t level = 4;
  size_t numAdditionalPoints = 0;
  bool normalizeWeights = false;

  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::CombiEvaluators::Collection evaluators(0);
  for (size_t d = 0; d < numDimensions; d++) {
    evaluators.push_back(sgpp::combigrid::CombiEvaluators::BSplineQuadrature(
        degree, weightFunctionsCollection[d], numAdditionalPoints, normalizeWeights));
  }
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::RegularLevelManager());
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);
  auto operation = std::make_shared<sgpp::combigrid::CombigridOperation>(
      pointHierarchies, evaluators, levelManager, gf, false);

  double integral = operation->evaluate(level);
  // int int (x^3+y^3)*(cos x*cos y)dxdy
  double exactSolution = 0.289025354482001;  // 1D: 0.1717381583560983;
  double error = fabs(integral - exactSolution);
  //  std::cout << "error: " << error << std::endl;
  BOOST_CHECK_SMALL(error, 1e-15);
}

double BSplineVarianceWithWeightsAndBounds(
    sgpp::combigrid::MultiIndex level, size_t degree, sgpp::combigrid::MultiFunction func,
    sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection,
    sgpp::base::DataVector bounds) {
  size_t numDimensions = level.size();
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());

  size_t numAdditionalPoints = 0;
  bool normalizeWeights = false;

  std::vector<
      std::shared_ptr<sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatArrayVector>>>
      evaluators(0);
  for (size_t d = 0; d < numDimensions; d++) {
    evaluators.push_back(std::make_shared<sgpp::combigrid::BSplineScalarProductEvaluator>(
        degree, weightFunctionsCollection[d], numAdditionalPoints, bounds[2 * d], bounds[2 * d + 1],
        normalizeWeights));
  }

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

// ToDo (rehmemk) What is wrong with this test? fix
// double oFunc(sgpp::base::DataVector const& v) { return std::pow(v[0], 3) + std::pow(v[1], 3); }
// double wFct(double x) { return sin(x); }
//
// BOOST_AUTO_TEST_CASE(testScalarProductsWithWeightFunctionAndBoundsOnLevel) {
//  // test on one level, for refinement
//  std::cout << "calculating mean and variance for f(x,y) = x^3+y^3 with weight function w(x) = "
//               "sin(x) on [0,2]^2 with B-splines of degree 3 "
//            << std::endl;
//  size_t numDimensions = 2;
//  size_t degree = 3;
//  sgpp::combigrid::MultiFunction func(oFunc);
//  sgpp::combigrid::SingleFunction weightfunction(wFct);
//  sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection(0);
//  sgpp::base::DataVector bounds;
//  for (size_t d = 0; d < numDimensions; d++) {
//    weightFunctionsCollection.push_back(weightfunction);
//    bounds.push_back(0);
//    bounds.push_back(2);
//  }
//
//  // ToDo (rehmemk) more levels like in the tests above
//  sgpp::combigrid::MultiIndex level{4, 4};
//  double variance =
//      BSplineVarianceWithWeightsAndBounds(level, degree, func, weightFunctionsCollection, bounds);
//
//  double realVariance = -0.575444693187592;
//  double varianceError = std::fabs(variance - realVariance);
//  //  std::cout << varianceError << std::endl;
//  BOOST_CHECK_SMALL(varianceError, 1e-13);
//}

#ifdef USE_DAKOTA

BOOST_AUTO_TEST_CASE(testBsplineStochasticCollocation_co2_lognormal) {
  std::cout << "Integrate objective function co2model and lognormal weight function  with B "
               "splines of degree 5 on level 4 in 1D "
            << std::endl;

  // create CO2 function and pdf weight functions
  sgpp::combigrid::CO2 co2Model;
  sgpp::combigrid::ProbabilityDensityFunction1DConfiguration pdf_config;
  pdf_config.pdfParameters.type_ =
      sgpp::combigrid::ProbabilityDensityFunctionType::BOUNDED_LOGNORMAL;
  pdf_config.pdfParameters.logmean_ = co2Model.logmean;
  pdf_config.pdfParameters.stddev_ = co2Model.stddev;
  pdf_config.pdfParameters.lowerBound_ = co2Model.bounds[0];
  pdf_config.pdfParameters.upperBound_ = co2Model.bounds[1];
  auto probabilityDensityFunction =
      std::make_shared<sgpp::combigrid::ProbabilityDensityFunction1D>(pdf_config);
  sgpp::combigrid::MultiFunction func(co2Model.eval);
  sgpp::combigrid::SingleFunction weight_function = probabilityDensityFunction->getWeightFunction();

  // initialize the Bspline surrogate model
  size_t numDims = co2Model.numDims;
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDims, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> storage;
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::AveragingLevelManager());
  sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection(0);
  for (size_t d = 0; d < numDims; d++) {
    weightFunctionsCollection.push_back(weight_function);
  }

  sgpp::combigrid::CombigridSurrogateModelConfiguration bsc_config;
  bsc_config.type = sgpp::combigrid::CombigridSurrogateModelsType::BSPLINE_STOCHASTIC_COLLOCATION;
  bsc_config.pointHierarchies = pointHierarchies;
  bsc_config.levelManager = levelManager;
  bsc_config.degree = 5;
  bsc_config.coefficientStorage = storage;
  bsc_config.weightFunctions = weightFunctionsCollection;
  bsc_config.bounds =
      sgpp::base::DataVector(std::vector<double>({co2Model.bounds[0], co2Model.bounds[1]}));

  sgpp::combigrid::BsplineStochasticCollocation bsc(bsc_config);

  // create level Structure and interpolate
  std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> newLevelStructure;
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> newCoefficientStorage;
  size_t numLevels = 4;
  sgpp::combigrid::GridFunction gf =
      BSplineCoefficientGridFunction(func, pointHierarchies, bsc_config.degree);
  createRegularLevelStructure(numLevels, bsc_config.degree, pointHierarchies, gf, false,
                              newLevelStructure, newCoefficientStorage, numDims);

  //  update config
  bsc_config.levelStructure = newLevelStructure;
  bsc_config.coefficientStorage = newCoefficientStorage;
  bsc.updateConfig(bsc_config);

  // check the moments
  //  std::cout << std::abs(co2Model.mean - bsc.mean()) << std::endl;
  //  std::cout << std::abs(co2Model.variance - bsc.variance()) << std::endl;
  BOOST_CHECK_SMALL(std::abs(co2Model.mean - bsc.mean()), 1e-8);
  BOOST_CHECK_SMALL(std::abs(co2Model.variance - bsc.variance()), 1e-9);
}

double objectiveFunction(sgpp::base::DataVector const& v) { return 4.0 * v[0] - 1.0; }

BOOST_AUTO_TEST_CASE(testBsplineNormalMeanAndVariance) {
  std::cout << "testing calculation of mean and variance using normally distribution and linear "
               "function transformed to [-1,3]"
            << std::endl;
  size_t numDimensions = 1;
  size_t degree = 5;
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(
      sgpp::combigrid::MultiFunction(objectiveFunction), pointHierarchies, degree);

  // set up the weight function collection as normally distributed probability density functions
  sgpp::combigrid::ProbabilityDensityFunction1DConfiguration pdf_config;
  pdf_config.pdfParameters.type_ = sgpp::combigrid::ProbabilityDensityFunctionType::NORMAL;
  pdf_config.pdfParameters.mean_ = 1.0;    // => we should obtain mean = 1.0
  pdf_config.pdfParameters.stddev_ = 0.1;  // => we should obtain variance = 0.01
  pdf_config.pdfParameters.lowerBound_ = -1;
  pdf_config.pdfParameters.upperBound_ = 3;
  auto probabilityDensityFunction =
      std::make_shared<sgpp::combigrid::ProbabilityDensityFunction1D>(pdf_config);
  sgpp::combigrid::SingleFunction oneDimensionsalWeightFunction =
      probabilityDensityFunction->getWeightFunction();
  sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection;
  sgpp::base::DataVector bounds;
  for (size_t d = 0; d < numDimensions; d++) {
    bounds.push_back(-1);
    bounds.push_back(3);
    weightFunctionsCollection.push_back(oneDimensionsalWeightFunction);
  }

  // set up the configuration for the B spline Stochastic Collocation
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> storage;
  sgpp::combigrid::CombigridSurrogateModelConfiguration config;
  config.type = sgpp::combigrid::CombigridSurrogateModelsType::BSPLINE_STOCHASTIC_COLLOCATION;
  config.pointHierarchies = pointHierarchies;
  config.levelManager = std::make_shared<sgpp::combigrid::AveragingLevelManager>();
  config.degree = degree;
  config.coefficientStorage = storage;
  config.weightFunctions = weightFunctionsCollection;
  config.bounds = bounds;
  sgpp::combigrid::BsplineStochasticCollocation bsc(config);

  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, config.degree);
  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));
  std::shared_ptr<sgpp::combigrid::LevelManager> varianceLevelManager(
      new sgpp::combigrid::RegularLevelManager());
  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  bool exploitNesting = false;
  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, Evaluators, varianceLevelManager, gf, exploitNesting,
      auxiliarySummationStrategyType);

  Operation->getLevelManager()->addRegularLevels(1);

  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager = Operation->getLevelManager();
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> coefficientStorage =
      Operation->getStorage();

  //  update config
  config.levelStructure = levelManager->getLevelStructure();
  config.coefficientStorage = coefficientStorage;
  bsc.updateConfig(config);

  double variance = bsc.variance();
  double ev = bsc.mean();
  //  std::cout << "mean: " << ev << " variance: " << variance << std::endl;
  BOOST_CHECK_SMALL(std::abs(ev - 1), 1e-13);
  BOOST_CHECK_SMALL(std::abs(variance - 0.01), 1e-13);
}

#endif

BOOST_AUTO_TEST_SUITE_END()
