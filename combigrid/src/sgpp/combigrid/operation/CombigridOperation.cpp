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
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridLinearCallbackEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PsiHermiteInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>

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
        fullGridEval(new FullGridLinearCallbackEvaluator<FloatScalarVector>(
            storage, evaluatorPrototypes, pointHierarchies)),
        combiEval(new CombigridEvaluator<FloatScalarVector>(pointHierarchies.size(), fullGridEval)),
        levelManager(levelManager) {
    levelManager->setLevelEvaluator(combiEval);
  }

  CombigridOperationImpl(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage,
      GridFunction gridFunc)
      : storage(storage),
        fullGridEval(new FullGridLinearGridBasedEvaluator<FloatScalarVector>(
            storage, evaluatorPrototypes, pointHierarchies, gridFunc)),
        combiEval(new CombigridEvaluator<FloatScalarVector>(pointHierarchies.size(), fullGridEval)),
        levelManager(levelManager) {
    levelManager->setLevelEvaluator(combiEval);
  }

  std::shared_ptr<AbstractCombigridStorage> storage;
  std::shared_ptr<AbstractFullGridEvaluator<FloatScalarVector>> fullGridEval;
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

  impl->fullGridEval->setParameters(scalars);
  impl->combiEval->clear();
}



std::shared_ptr<AbstractFullGridEvaluator<FloatScalarVector>>  CombigridOperation::getFullGridEvaluator(){
return impl->fullGridEval;

}


double CombigridOperation::getResult() { return impl->combiEval->getValue().value(); }

double CombigridOperation::evaluate(size_t q, base::DataVector const& param) {
  setParameters(param);
  impl->levelManager->addRegularLevels(q);

  return getResult();
}

std::shared_ptr<AbstractCombigridStorage> CombigridOperation::getStorage() { return impl->storage; }

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

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformBoundaryPsiLinearInterpolation(
    size_t numDimensions, size_t psiDimension, MultiFunction func) {
  std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluators(
      numDimensions, CombiEvaluators::linearInterpolation());

  evaluators[psiDimension] = CombiEvaluators::psiHermiteInterpolation();

  return std::make_shared<CombigridOperation>(std::vector<std::shared_ptr<AbstractPointHierarchy>>(
                                                  numDimensions, CombiHierarchies::expUniformBoundary()),
                                              evaluators, std::make_shared<StandardLevelManager>(),
                                              func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformBoundaryZetaLinearInterpolation(
    size_t numDimensions, size_t zetaDimension, MultiFunction func) {
  std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluators(
      numDimensions, CombiEvaluators::linearInterpolation());

  evaluators[zetaDimension] = CombiEvaluators::zetaHermiteInterpolation();

  return std::make_shared<CombigridOperation>(std::vector<std::shared_ptr<AbstractPointHierarchy>>(
                                                  numDimensions, CombiHierarchies::expUniformBoundary()),
                                              evaluators, std::make_shared<StandardLevelManager>(),
                                              func);
}




std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformPsiHermiteInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniform()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::psiHermiteInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformBoundaryPsiHermiteInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniformBoundary()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::psiHermiteInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformPsiHermiteInterpolation(
    size_t numDimensions, size_t zetaDimension, MultiFunction func) {
  std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluators(
      numDimensions, CombiEvaluators::psiHermiteInterpolation());

  evaluators[zetaDimension] = CombiEvaluators::zetaHermiteInterpolation();

  return std::make_shared<CombigridOperation>(std::vector<std::shared_ptr<AbstractPointHierarchy>>(
                                                  numDimensions, CombiHierarchies::expUniform()),
                                              evaluators, std::make_shared<StandardLevelManager>(),
                                              func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformBoundaryPsiHermiteInterpolation(
    size_t numDimensions, size_t zetaDimension, MultiFunction func) {
  std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluators(
      numDimensions, CombiEvaluators::psiHermiteInterpolation());

  evaluators[zetaDimension] = CombiEvaluators::zetaHermiteInterpolation();

  return std::make_shared<CombigridOperation>(std::vector<std::shared_ptr<AbstractPointHierarchy>>(
                                                  numDimensions, CombiHierarchies::expUniformBoundary()),
                                              evaluators, std::make_shared<StandardLevelManager>(),
                                              func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformBoundaryZetaHermiteInterpolation(
    size_t numDimensions, size_t zetaDimension, MultiFunction func) {
  std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluators(
      numDimensions, CombiEvaluators::psiHermiteInterpolation());

  evaluators[zetaDimension] = CombiEvaluators::zetaHermiteInterpolation();

  return std::make_shared<CombigridOperation>(std::vector<std::shared_ptr<AbstractPointHierarchy>>(
                                                  numDimensions, CombiHierarchies::expUniformBoundary()),
                                              evaluators, std::make_shared<StandardLevelManager>(),
                                              func);
}


std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformZetaHermiteInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniform()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::zetaHermiteInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}


std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformBoundaryZetaHermiteInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniformBoundary()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>(
          numDimensions, CombiEvaluators::zetaHermiteInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}


// TODO andere psi funktion anpassen
std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformBoundaryZetaInterpolation(
    size_t numDimensions, std::vector<int> zetaDimensions, MultiFunction func) {
  std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluators(
      numDimensions, CombiEvaluators::psiHermiteInterpolation());

      for(int i : zetaDimensions) {
        evaluators[i]=CombiEvaluators::zetaHermiteInterpolation();
      }


  return std::make_shared<CombigridOperation>(std::vector<std::shared_ptr<AbstractPointHierarchy>>(
                                                  numDimensions, CombiHierarchies::expUniformBoundary()),
                                              evaluators, std::make_shared<StandardLevelManager>(),
                                              func);
}

std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformBoundaryPsiInterpolation(
    size_t numDimensions, std::vector<int> psiDimensions, MultiFunction func) {
  std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluators(
      numDimensions, CombiEvaluators::zetaHermiteInterpolation());

      for(int i : psiDimensions) {
        evaluators[i]=CombiEvaluators::psiHermiteInterpolation();
      }


  return std::make_shared<CombigridOperation>(std::vector<std::shared_ptr<AbstractPointHierarchy>>(
                                                  numDimensions, CombiHierarchies::expUniformBoundary()),
                                              evaluators, std::make_shared<StandardLevelManager>(),
                                              func);
}


// ToDo (rehmemk) Choose basistype SLE here! instead of first linear and changing to SLE in
// AbstractFullGridinearEvaluator.hpp
std::shared_ptr<CombigridOperation> CombigridOperation::createExpUniformBsplineInterpolation(
    size_t numDimensions, MultiFunction func) {
  return std::make_shared<CombigridOperation>(
      std::vector<std::shared_ptr<AbstractPointHierarchy>>(numDimensions,
                                                           CombiHierarchies::expUniform()),
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>>( numDimensions, CombiEvaluators::linearInterpolation()),
      std::make_shared<StandardLevelManager>(), func);
}

 

} /* namespace combigrid */
} /* namespace sgpp*/
