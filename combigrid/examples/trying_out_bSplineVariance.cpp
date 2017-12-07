// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/common/GridConversion.hpp>
#include <sgpp/combigrid/integration/MCIntegrator.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/OperationConfiguration.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/RegularLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridGridBasedEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridLinearSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridQuadraticSummationStrategy.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineQuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineScalarProductEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotNakBsplineBoundaryCombigrid.hpp>
#include <sgpp/quadrature/sampling/NaiveSampleGenerator.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

double atanMean = 3.514491266445763;
double atanMeanSquare = 15.804752455860012;
double atanVariance = 3.453103593936470;

double xcubeMean = 0.25;
double xcubeMeanSquare = 0.142857142857143;
double xcubeVariance = 0.080357142857143;

double xcube2Mean = 0.500000000000164;
double xcube2MeanSquare = 0.410714285751263;
double xcube2Variance = 0.160714285751100;

double sinexpMean = 0.672377634914871;
double sinexpMeanSquare = 0.644728986265339;
double sinexpVariance = 0.192637302331624;

double xquintMean = 1.0 / 6.0;
double xquintMeanSquare = 1.0 / 11.0;
double xquintVariance = 1.0 / 11.0 - 1.0 / 36.0;

size_t numDimensions = 2;
double f(sgpp::base::DataVector const& v) {
  //  return v[0] * v[0] * v[0] * v[0] * v[0];
  //  return v[0] * v[0] * v[0] + v[1] * v[1] * v[1];
  //  return std::sin(v[0]) * std::exp(v[1] * v[1]);
  return std::atan(50 * (v[0] - .35)) + M_PI / 2 + 4 * std::pow(v[1], 3) +
         std::exp(v[0] * v[1] - 1);
}
void printLevelstructure(
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelstructure) {
  auto it = levelstructure->getStoredDataIterator();
  while (it->isValid()) {
    sgpp::combigrid::MultiIndex index = it->getMultiIndex();
    for (auto& i : index) {
      std::cout << i << " ";
    }
    std::cout << "\n";
    it->moveToNext();
  }
}

void printLevelstructurePoints(
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelStructure,
    sgpp::combigrid::CombiHierarchies::Collection pointHierarchies,
    std::vector<bool> orderingConfiguration) {
  auto levelIterator = levelStructure->getStoredDataIterator();
  while (levelIterator->isValid()) {
    sgpp::combigrid::MultiIndex level = levelIterator->getMultiIndex();
    for (size_t d = 0; d < numDimensions; ++d) {
      std::cout << "level " << level[d] << std::endl;
      auto gridpoints = pointHierarchies[d]->getPoints(level[d], orderingConfiguration[d]);
      std::cout << "point hierarchies: ";
      for (auto& p : gridpoints) {
        std::cout << p << " ";
      }
      std::cout << "\n";
    }
    levelIterator->moveToNext();
  }
}

void printSGGridToFile(sgpp::base::GridStorage const& gridStorage) {
  std::string plotstr = "/home/rehmemk/SGS_Sync/Plotting/combigrid_bsplines/convertedGrid.dat";
  remove(plotstr.c_str());
  std::ofstream plotfile;
  plotfile.open(plotstr.c_str(), std::ios::app);
  plotfile << "#grid points" << std::endl;
  for (size_t q = 0; q < gridStorage.getSize(); q++) {
    auto point = gridStorage.getPoint(q);
    for (size_t d = 0; d < numDimensions - 1; d++) {
      plotfile << point.getStandardCoordinate(d) << ", ";
    }
    plotfile << point.getStandardCoordinate(numDimensions - 1);
    plotfile << "\n";
  }
  plotfile.close();
}

void calculateCGerror(
    double& maxErr, double& L2Err,
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> Operation,
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelStructure) {
  maxErr = 0;
  L2Err = 0;
  size_t numMCpoints = 10000;
  sgpp::quadrature::NaiveSampleGenerator generator(numDimensions);
  sgpp::base::DataMatrix params(numDimensions, numMCpoints);
  sgpp::base::DataVector p(numDimensions);
  sgpp::base::DataVector Feval(numMCpoints, 0.0);
  for (size_t i = 0; i < numMCpoints; i++) {
    generator.getSample(p);
    params.setColumn(i, p);
    Feval.set(i, f(p));
  }
  Operation->setParameters(params);
  Operation->getLevelManager()->addLevelsFromStructure(levelStructure);
  sgpp::base::DataVector CGeval = Operation->getResult();
  CGeval.sub(Feval);
  for (size_t i = 0; i < CGeval.size(); i++) {
    maxErr = (fabs(CGeval[i]) > maxErr) ? fabs(CGeval[i]) : maxErr;
    L2Err += fabs(CGeval[i] * CGeval[i]);
  }
  L2Err = sqrt(L2Err / static_cast<double>(numMCpoints));
}

void calculateSGerror(double& maxErr, double& L2Err,
                      sgpp::optimization::InterpolantScalarFunction& u) {
  double diff = 0;
  maxErr = 0;
  L2Err = 0;
  size_t numMCpoints = 10000;
  sgpp::quadrature::NaiveSampleGenerator generator(numDimensions);
  sgpp::base::DataVector p(numDimensions);
  for (size_t i = 0; i < numMCpoints; i++) {
    generator.getSample(p);
    diff = fabs(u.eval(p) - f(p));
    maxErr = (diff > maxErr) ? diff : maxErr;
    L2Err += diff * diff;
  }
  L2Err = sqrt(L2Err / static_cast<double>(numMCpoints));
}

void calculateCGSGDifference(
    double& maxErr, double& L2Err,
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> Operation,
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelStructure,
    sgpp::optimization::InterpolantScalarFunction& u) {
  maxErr = 0;
  L2Err = 0;
  size_t numMCpoints = 10000;
  sgpp::quadrature::NaiveSampleGenerator generator(numDimensions);
  sgpp::base::DataMatrix params(numDimensions, numMCpoints);
  sgpp::base::DataVector p(numDimensions);
  sgpp::base::DataVector SGeval(numMCpoints, 0.0);
  for (size_t i = 0; i < numMCpoints; i++) {
    generator.getSample(p);
    params.setColumn(i, p);
    SGeval.set(i, u.eval(p));
  }
  Operation->setParameters(params);
  Operation->getLevelManager()->addLevelsFromStructure(levelStructure);
  sgpp::base::DataVector CGeval = Operation->getResult();
  CGeval.sub(SGeval);
  for (size_t i = 0; i < CGeval.size(); i++) {
    maxErr = (fabs(CGeval[i]) > maxErr) ? fabs(CGeval[i]) : maxErr;
    L2Err += fabs(CGeval[i] * CGeval[i]);
  }
  L2Err = sqrt(L2Err / static_cast<double>(numMCpoints));
}

std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> createVarianceLevelStructure(
    size_t const& numlevels, size_t const& degree,
    sgpp::combigrid::CombiHierarchies::Collection const& pointHierarchies,
    sgpp::combigrid::GridFunction gf, bool exploitNesting) {
  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineScalarProduct, degree);
  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::AveragingLevelManager());
  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::VARIANCE;

  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, Evaluators, levelManager, gf, exploitNesting,
      auxiliarySummationStrategyType);

  // create level structure
  // First step is to guarantee the existence of level (1,..,1). Otherwise a conversion to an SG
  // grid wouldn't be possible
  Operation->getLevelManager()->addRegularLevels(numDimensions + numlevels);
  // ToDo(rehmemk) the number of levels added is now a little mysterious because the user does not
  // know how many levels are added to guarantee (1,..,1) before the number of levels he chose to
  // be added is added
  //  Operation->getLevelManager()->addLevelsAdaptiveByNumLevels(numlevels);
  auto levelStructure = Operation->getLevelManager()->getLevelStructure();
  return levelStructure;
}

void BSplineGridConversion(size_t degree, size_t numlevels) {
  // create interpolation operation
  sgpp::combigrid::MultiFunction func(f);
  sgpp::combigrid::EvaluatorConfiguration evalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, degree);
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::CombiEvaluators::MultiCollection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(evalConfig));
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  bool exploitNesting = false;
  std::shared_ptr<sgpp::combigrid::LevelManager> dummyLevelManager(
      new sgpp::combigrid::RegularLevelManager());
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);
  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, evaluators, dummyLevelManager, gf, exploitNesting, summationStrategyType);

  // create variance adaptive level structure
  auto levelStructure =
      createVarianceLevelStructure(numlevels, degree, pointHierarchies, gf, exploitNesting);

  std::vector<bool> orderingConfiguration;
  for (size_t d = 0; d < numDimensions; ++d) {
    orderingConfiguration.push_back(evaluators[d]->needsOrderedPoints());
  }

  // convert level structure to SG
  std::shared_ptr<sgpp::base::Grid> grid;
  grid.reset(sgpp::base::Grid::createNakBsplineBoundaryCombigridGrid(numDimensions, degree));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  convertexpUniformBoundaryCombigridToHierarchicalSparseGrid(levelStructure, gridStorage);

  //  print options
  //  printLevelstructure(levelStructure);
  //  printSGGridToFile(gridStorage);

  // interpolate on SG
  sgpp::base::DataMatrix interpolParams(numDimensions, gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
      p[j] = gp.getStandardCoordinate(j);
    }
    interpolParams.setColumn(i, p);
  }

  // obtain function values from combigrid surrogate
  Operation->setParameters(interpolParams);
  Operation->getLevelManager()->addLevelsFromStructure(levelStructure);
  sgpp::base::DataVector f_values = Operation->getResult();

  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::optimization::sle_solver::Auto sleSolver;
  sgpp::base::DataVector alpha(grid->getSize());
  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
    std::cout << "Solving failed!" << std::endl;
  }
  sgpp::optimization::InterpolantScalarFunction u(*grid, alpha);

  //  std::cout << "num CG points: " << Operation->getLevelManager()->numGridPoints();
  //  std::cout << ", num SG points " << gridStorage.getSize() << std::endl;
  std::cout << gridStorage.getSize() << " ";

  // error calculations
  double CGL2Err, CGMaxErr, SGL2Err, SGMaxErr, CompL2Err, CompMaxErr = 0;
  calculateCGerror(CGMaxErr, CGL2Err, Operation, levelStructure);
  calculateSGerror(SGMaxErr, SGL2Err, u);
  calculateCGSGDifference(CompMaxErr, CompL2Err, Operation, levelStructure, u);

  //  std::cout << "\n";
  //  std::cout << "CG L2:   " << CGL2Err << "   CG max:   " << CGMaxErr << std::endl;
  //  std::cout << "SG L2:   " << SGL2Err << "   SG max:   " << SGMaxErr << std::endl;
  //  std::cout << "Comp L2: " << CompL2Err << " Comp max: " << CompMaxErr << std::endl;
  //  std::cout << "L2 error: " << CGL2Err << " ";

  sgpp::base::DataVector dummyparams;
  auto quadOperation =
      sgpp::combigrid::CombigridOperation::createExpUniformBoundaryBsplineQuadrature(numDimensions,
                                                                                     func, degree);
  quadOperation->getLevelManager()->addLevelsFromStructure(levelStructure);
  double mean = quadOperation->getResult();

  sgpp::base::Grid* gridptr = grid.get();
  sgpp::pde::OperationMatrixLTwoDotNakBsplineBoundaryCombigrid massMatrix(gridptr);
  sgpp::base::DataVector product(alpha.size(), 0);
  massMatrix.mult(alpha, product);
  double meanSquare = product.dotProduct(alpha);
  double variance = meanSquare - mean * mean;
  //  std::cout << " mean error: " << fabs(mean - atanMean) << " ";
  //  std::cout << " meanSquare error : " << fabs(meanSquare - atanMeanSquare) << " ";
  std::cout << fabs(variance - atanVariance) << std::endl;

  // MATRIX DEBUGGING:

  //  //  for (size_t i = 0; i < alpha.getSize(); i++) {
  //  //    for (size_t j = 0; j < alpha.getSize(); j++) {
  //  sgpp::base::DataVector ivec(alpha.getSize(), 0.0);
  //  sgpp::base::DataVector jvec(alpha.getSize(), 0.0);
  //  //      ivec[i] = 1;
  //  //      jvec[j] = 1;
  //  ivec[0] = 1;
  //  jvec[0] = 1;
  //  massMatrix.mult(ivec, product);
  //  std::cout << std::fixed << product.dotProduct(jvec) << " ";
  //  //}
  //  std::cout << "\n";
  //  //  }
}

int main() {
  size_t degree = 5;
  size_t numLevels = 8;
  for (size_t maxLevel = 0; maxLevel < numLevels; maxLevel++) {
    //    std::cout << maxLevel << ", ";
    //  size_t maxLevel = 0;
    BSplineGridConversion(degree, maxLevel);
  }
  return 0;
}
