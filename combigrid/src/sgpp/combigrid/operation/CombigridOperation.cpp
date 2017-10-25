// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/UniformPointDistribution.hpp>
#include <sgpp/combigrid/grid/growth/ExponentialGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridLinearSummationStrategy.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <vector>

typedef sgpp::combigrid::AveragingLevelManager StandardLevelManager;

namespace sgpp {
namespace combigrid {

class CombigridOperationImpl {
 public:
  CombigridOperationImpl(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage)
      : storage(storage),
        pointHierarchies(pointHierarchies),
        summationStrategy(new FullGridLinearSummationStrategy<FloatScalarVector>(
            storage, evaluatorPrototypes, pointHierarchies)),
        combiEval(
            new CombigridEvaluator<FloatScalarVector>(pointHierarchies.size(), summationStrategy)),

        levelManager(levelManager) {
    levelManager->setLevelEvaluator(combiEval);
  }

  CombigridOperationImpl(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage,
      GridFunction gridFunc)
      : storage(storage),
        pointHierarchies(pointHierarchies),
        summationStrategy(new FullGridLinearSummationStrategy<FloatScalarVector>(
            storage, evaluatorPrototypes, pointHierarchies, gridFunc)),
        combiEval(
            new CombigridEvaluator<FloatScalarVector>(pointHierarchies.size(), summationStrategy)),
        levelManager(levelManager) {
    levelManager->setLevelEvaluator(combiEval);
  }

  std::shared_ptr<AbstractCombigridStorage> storage;
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies;
  std::shared_ptr<AbstractFullGridSummationStrategy<FloatScalarVector>> summationStrategy;
  std::shared_ptr<CombigridEvaluator<FloatScalarVector>> combiEval;
  std::shared_ptr<LevelManager> levelManager;
};

CombigridOperation::CombigridOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
    std::shared_ptr<LevelManager> levelManager, MultiFunction func)
    : impl(new CombigridOperationImpl(pointHierarchies, evaluatorPrototypes, levelManager,
                                      std::shared_ptr<AbstractCombigridStorage>(
                                          new CombigridTreeStorage(pointHierarchies, func)))) {}

CombigridOperation::CombigridOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
    std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage)
    : impl(new CombigridOperationImpl(pointHierarchies, evaluatorPrototypes, levelManager,
                                      storage)) {}

CombigridOperation::CombigridOperation(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
    std::shared_ptr<LevelManager> levelManager, GridFunction gridFunc, bool exploitNesting)
    : impl(new CombigridOperationImpl(
          pointHierarchies, evaluatorPrototypes, levelManager,
          std::shared_ptr<AbstractCombigridStorage>(
              new CombigridTreeStorage(pointHierarchies, exploitNesting)),
          gridFunc)) {}

void CombigridOperation::setParameters(const base::DataVector& param) {
  std::vector<FloatScalarVector> scalars(param.getSize());
  for (size_t i = 0; i < param.getSize(); ++i) {
    scalars[i].value() = param[i];
  }

  impl->summationStrategy->setParameters(scalars);
  impl->combiEval->clear();
}

double CombigridOperation::getResult() { return impl->combiEval->getValue().value(); }

double CombigridOperation::evaluate(size_t q, base::DataVector const& param) {
  setParameters(param);
  impl->levelManager->addRegularLevels(q);

  return getResult();
}

std::shared_ptr<AbstractCombigridStorage> CombigridOperation::getStorage() { return impl->storage; }
std::vector<std::shared_ptr<AbstractPointHierarchy>> CombigridOperation::getPointHierarchies() {
  return impl->pointHierarchies;
}

std::shared_ptr<LevelManager> CombigridOperation::getLevelManager() { return impl->levelManager; }

void CombigridOperation::setLevelManager(std::shared_ptr<LevelManager> levelManager) {
  levelManager->setLevelEvaluator(impl->combiEval);
  impl->levelManager = levelManager;
}

size_t CombigridOperation::numStoredFunctionValues() { return impl->storage->getNumEntries(); }

size_t CombigridOperation::numGridPoints() { return impl->levelManager->numGridPoints(); }

size_t CombigridOperation::getUpperPointBound() const {
  return impl->levelManager->getUpperPointBound();
}

std::shared_ptr<CombigridOperation>
CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(size_t numDimensions,
                                                                   MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expClenshawCurtis()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpChebyshevPolynomialInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expChebyshev()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpLejaPolynomialInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expLeja()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpL2LejaPolynomialInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expL2Leja()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformPolynomialInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniform()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation>
CombigridOperation::createExpUniformBoundaryPolynomialInterpolation(size_t numDimensions,
                                                                    MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniformBoundary()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformLinearInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniform()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::linearInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformBoundaryLinearInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniformBoundary()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::linearInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation>
CombigridOperation::createLinearClenshawCurtisPolynomialInterpolation(size_t numDimensions,
                                                                      MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, std::make_shared<NonNestedPointHierarchy>(
                             std::make_shared<ClenshawCurtisDistribution>(),
                             std::make_shared<IdentityPointOrdering>(
                                 std::make_shared<LinearGrowthStrategy>(2), true))),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearLejaPolynomialInterpolation(
    size_t numDimensions, MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearLeja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearL2LejaPolynomialInterpolation(
    size_t numDimensions, MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearL2Leja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearUniformPolynomialInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, std::make_shared<NonNestedPointHierarchy>(
                             std::make_shared<UniformPointDistribution>(),
                             std::make_shared<IdentityPointOrdering>(
                                 std::make_shared<LinearGrowthStrategy>(2), true))),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation>
CombigridOperation::createLinearUniformBoundaryPolynomialInterpolation(size_t numDimensions,
                                                                       MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, std::make_shared<NonNestedPointHierarchy>(
                             std::make_shared<UniformBoundaryPointDistribution>(),
                             std::make_shared<IdentityPointOrdering>(
                                 std::make_shared<LinearGrowthStrategy>(2), true))),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::polynomialInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearLejaQuadrature(
    size_t numDimensions, MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearLeja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::quadrature()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearL2LejaQuadrature(
    size_t numDimensions, MultiFunction func, size_t growthFactor) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(
          numDimensions, CombiHierarchies::linearL2Leja(growthFactor)),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::quadrature()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpClenshawCurtisQuadrature(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expClenshawCurtis()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::quadrature()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> auxiliaryBsplineFunction(size_t numDimensions,
                                                             MultiFunction func,
                                                             size_t hierarchyType,
                                                             size_t operationType,
                                                             size_t growthFactor, size_t degree) {
  // ToDo (rehmemk) Ersetze hierarchyType und operationType durch Übergabe der entsprechenden
  // Collections, dies ermöglicht insbesondere gemischte Evaluatoren
  sgpp::combigrid::CombiHierarchies::Collection grids;
  if (hierarchyType == 1) {
    grids.resize(numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  } else if (hierarchyType == 2) {
    grids.resize(numDimensions, sgpp::combigrid::CombiHierarchies::expClenshawCurtis());
  } else if (hierarchyType == 3) {
    grids.resize(numDimensions, sgpp::combigrid::CombiHierarchies::expChebyshev());
  } else if (hierarchyType == 4) {
    grids.resize(numDimensions, sgpp::combigrid::CombiHierarchies::linearLeja(growthFactor));
  } else if (hierarchyType == 5) {
    grids.resize(numDimensions, sgpp::combigrid::CombiHierarchies::linearL2Leja(growthFactor));
  }

  sgpp::combigrid::CombiEvaluators::Collection evaluators;
  if (operationType == 1) {
    evaluators.resize(numDimensions,
                      sgpp::combigrid::CombiEvaluators::BSplineInterpolation(degree));
  } else if (operationType == 2) {
    evaluators.resize(numDimensions, sgpp::combigrid::CombiEvaluators::BSplineQuadrature(degree));
  }
  // So far only WeightedRatioLevelManager has been used
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::WeightedRatioLevelManager());

  // stores the values of the objective function
  auto funcStorage = std::make_shared<sgpp::combigrid::CombigridTreeStorage>(grids, func);

  // Grid Function that calculates the coefficients for the B-Spline interpolation.
  // The coeficients for each B-Spline are saved in a TreeStorage encoded by a MultiIndex
  sgpp::combigrid::GridFunction gf([=](std::shared_ptr<sgpp::combigrid::TensorGrid> grid) {
    sgpp::combigrid::CombiEvaluators::Collection interpolEvaluators(
        numDimensions, sgpp::combigrid::CombiEvaluators::BSplineInterpolation(degree));
    size_t numDimensions = grid->getDimension();
    auto coefficientTree = std::make_shared<sgpp::combigrid::TreeStorage<double>>(numDimensions);
    auto level = grid->getLevel();
    std::vector<size_t> numGridPointsVec = grid->numPoints();
    size_t numGridPoints = 1;
    for (size_t i = 0; i < numGridPointsVec.size(); i++) {
      numGridPoints *= numGridPointsVec[i];
    }

    sgpp::combigrid::CombiEvaluators::Collection evalCopy(numDimensions);
    for (size_t dim = 0; dim < numDimensions; ++dim) {
      evalCopy[dim] = interpolEvaluators[dim]->cloneLinear();
      bool needsSorted = evalCopy[dim]->needsOrderedPoints();
      auto gridPoints = grids[dim]->getPoints(level[dim], needsSorted);
      evalCopy[dim]->setGridPoints(gridPoints);
    }
    sgpp::base::DataMatrix A(numGridPoints, numGridPoints);
    sgpp::base::DataVector coefficients_sle(numGridPoints);
    sgpp::base::DataVector functionValues(numGridPoints);

    // Creates an iterator that yields the multi-indices of all grid points in the grid.
    sgpp::combigrid::MultiIndexIterator it(grid->numPoints());
    auto funcIter =
        funcStorage->getGuidedIterator(level, it, std::vector<bool>(numDimensions, true));

    for (size_t ixEvalPoints = 0; funcIter->isValid(); ++ixEvalPoints, funcIter->moveToNext()) {
      auto gridPoint = grid->getGridPoint(funcIter->getMultiIndex());
      functionValues[ixEvalPoints] = funcIter->value();

      std::vector<std::vector<double>> basisValues;
      for (size_t dim = 0; dim < numDimensions; ++dim) {
        evalCopy[dim]->setParameter(sgpp::combigrid::FloatScalarVector(gridPoint[dim]));
        auto basisValues1D = evalCopy[dim]->getBasisValues();
        // basis values at gridPoint
        std::vector<double> basisValues1D_vec(basisValues1D.size());
        for (size_t i = 0; i < basisValues1D.size(); i++) {
          basisValues1D_vec[i] = basisValues1D[i].value();
        }
        basisValues.push_back(basisValues1D_vec);
      }

      sgpp::combigrid::MultiIndexIterator innerIter(grid->numPoints());
      for (size_t ixBasisFunctions = 0; innerIter.isValid();
           ++ixBasisFunctions, innerIter.moveToNext()) {
        double splineValue = 1.0;
        auto innerIndex = innerIter.getMultiIndex();
        for (size_t dim = 0; dim < numDimensions; ++dim) {
          splineValue *= basisValues[dim][innerIndex[dim]];
        }
        A.set(ixEvalPoints, ixBasisFunctions, splineValue);
      }
    }

    sgpp::optimization::FullSLE sle(A);
    sgpp::optimization::sle_solver::Auto solver;
    sgpp::optimization::Printer::getInstance().setVerbosity(-1);
    bool solved = solver.solve(sle, functionValues, coefficients_sle);

    /*std::cout << A.toString() << std::endl;
    std::cout << "fct: ";
    for (size_t i = 0; i < functionValues.size(); i++) {
      std::cout << functionValues[i] << " ";
    }
    std::cout << "\ncoeff: ";
    for (size_t i = 0; i < coefficients_sle.size(); i++) {
      std::cout << coefficients_sle[i] << " ";
    }
    std::cout << "\n";
    std::cout << "--------" << std::endl;
    */
    if (!solved) {
      exit(-1);
    }

    it.reset();
    for (size_t vecIndex = 0; it.isValid(); ++vecIndex, it.moveToNext()) {
      coefficientTree->set(it.getMultiIndex(), coefficients_sle[vecIndex]);
    }

    return coefficientTree;
  });

  /**
   * We have to specify if the function always produces the same value for the same grid points.
   * This can make the storage smaller if the grid points are nested. In this implementation,
   * this is false.
   */
  bool exploitNesting = false;
  return std::make_shared<sgpp::combigrid::CombigridOperation>(grids, evaluators, levelManager, gf,
                                                               exploitNesting);
}

std::shared_ptr<CombigridOperation>
CombigridOperation::createExpUniformBoundaryBsplineInterpolation(size_t numDimensions,
                                                                 MultiFunction func,
                                                                 size_t degree = 3) {
  size_t dummygrowthfactor = 0;
  size_t gridType = 1;
  size_t operationType = 1;
  return auxiliaryBsplineFunction(numDimensions, func, gridType, operationType, dummygrowthfactor,
                                  degree);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpClenshawCurtisBsplineInterpolation(
    size_t numDimensions, MultiFunction func, size_t degree = 3) {
  size_t dummygrowthfactor = 0;
  size_t gridType = 2;
  size_t operationType = 1;
  return auxiliaryBsplineFunction(numDimensions, func, gridType, operationType, dummygrowthfactor,
                                  degree);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpChebyshevBsplineInterpolation(
    size_t numDimensions, MultiFunction func, size_t degree = 3) {
  size_t dummygrowthfactor = 0;
  size_t gridType = 3;
  size_t operationType = 1;
  return auxiliaryBsplineFunction(numDimensions, func, gridType, operationType, dummygrowthfactor,
                                  degree);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearLejaBsplineInterpolation(
    size_t numDimensions, MultiFunction func, size_t degree = 3, size_t growthFactor = 2) {
  size_t gridType = 4;
  size_t operationType = 1;
  return auxiliaryBsplineFunction(numDimensions, func, gridType, operationType, growthFactor,
                                  degree);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createLinearL2LejaBsplineInterpolation(
    size_t numDimensions, MultiFunction func, size_t degree = 3, size_t growthFactor = 2) {
  size_t gridType = 5;
  size_t operationType = 1;
  return auxiliaryBsplineFunction(numDimensions, func, gridType, operationType, growthFactor,
                                  degree);
}
std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformBoundaryBsplineQuadrature(
    size_t numDimensions, MultiFunction func, size_t degree = 3) {
  size_t dummygrowthfactor = 0;
  size_t gridType = 1;
  size_t operationType = 2;
  return auxiliaryBsplineFunction(numDimensions, func, gridType, operationType, dummygrowthfactor,
                                  degree);
}

} /* namespace combigrid */
} /* namespace sgpp*/
