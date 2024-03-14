// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingRelativeDensityRatioEstimation.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/algorithm/SystemMatrixRelativeDensityRatioEstimation.hpp>
#include <sgpp/solver/SLESolver.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusAbsValueCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusAbsValueRefinementFunctor.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

#include <string>
#include <vector>

using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::CoarseningFunctor;
using sgpp::base::CoarseningFunctorType;
using sgpp::base::RefinementFunctor;
using sgpp::base::RefinementFunctorType;
using sgpp::base::SurplusCoarseningFunctor;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::base::SurplusVolumeCoarseningFunctor;
using sgpp::base::SurplusVolumeRefinementFunctor;
using sgpp::base::SurplusAbsValueCoarseningFunctor;
using sgpp::base::SurplusAbsValueRefinementFunctor;

using sgpp::base::application_exception;

using sgpp::solver::SLESolver;

namespace sgpp {
namespace datadriven {

ModelFittingRelativeDensityRatioEstimation::ModelFittingRelativeDensityRatioEstimation(
    const FitterConfigurationDensityLeastSquares &config)
    : ModelFittingBaseSingleGrid{}, refinementsPerformed{0}, initialGridSize{0} {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityLeastSquares>(config));
  solver = std::unique_ptr<SLESolver>{buildSolver(this->config->getSolverFinalConfig())};
  std::cout << "omega: " << this->config->getDensityEstimationConfig().omega_ << std::endl;
}

// TODO(lettrich): exceptions have to be thrown if not valid.
double ModelFittingRelativeDensityRatioEstimation::evaluate(const DataVector &sample) {
  auto opEval = std::unique_ptr<base::OperationEval>{op_factory::createOperationEval(*grid)};
  return opEval->eval(alpha, sample);
}

// TODO(lettrich): exceptions have to be thrown if not valid.
void ModelFittingRelativeDensityRatioEstimation::evaluate(DataMatrix &samples,
                                                          DataVector &results) {
  auto opMultEval = std::unique_ptr<base::OperationMultipleEval>{
      op_factory::createOperationMultipleEval(*grid, samples, config->getMultipleEvalConfig())};
  opMultEval->eval(alpha, results);
}

double ModelFittingRelativeDensityRatioEstimation::relativeLeastSquaresLossApprox(
    DataMatrix &samplesP, DataMatrix &samplesQ) {
  double numP = static_cast<double>(samplesP.getNrows());
  double numQ = static_cast<double>(samplesQ.getNrows());
  DataVector fp(samplesP.getNrows());
  DataVector fq(samplesQ.getNrows());
  this->evaluate(samplesP, fp);
  this->evaluate(samplesQ, fq);
  double om = config->getDensityEstimationConfig().omega_;  // just to make it cleaner
  return om * fp.sumsqr() / numP + (1. - om) * fq.sumsqr() / numQ - 2 * fp.sum() / numP;
}

double ModelFittingRelativeDensityRatioEstimation::relativeKLDivergenceApprox(
    DataMatrix &samplesP) {
  // samplesQ not needed to compute this
  double numP = static_cast<double>(samplesP.getNrows());
  DataVector fp(samplesP.getNrows());
  this->evaluate(samplesP, fp);
  return fp.sumlog() / numP;
}

double ModelFittingRelativeDensityRatioEstimation::relativePEDivergenceApprox(
    DataMatrix &samplesP, DataMatrix &samplesQ) {
  double numP = static_cast<double>(samplesP.getNrows());
  double numQ = static_cast<double>(samplesQ.getNrows());
  DataVector fp(samplesP.getNrows());
  DataVector fq(samplesQ.getNrows());
  this->evaluate(samplesP, fp);
  this->evaluate(samplesQ, fq);
  double om = config->getDensityEstimationConfig().omega_;  // just to make it cleaner
  return -om / 2. * fp.sumsqr() / numP - (1. - om) / 2. * fq.sumsqr() / numQ + fp.sum() / numP -
         1. / 2.;
}

void ModelFittingRelativeDensityRatioEstimation::fit(Dataset &newDatasetP, Dataset &newDatasetQ) {
  // clear model
  reset();
  this->dataset = &newDatasetP;
  this->extraDataset = &newDatasetQ;

  // build grid
  auto &gridConfig = config->getGridConfig();
  gridConfig.dim_ = dataset->getDimension();  // extraDataset works as well
  grid = std::unique_ptr<Grid>{buildGrid(config->getGridConfig())};
  // build surplus vector
  alpha = DataVector(grid->getSize());

  assembleSystemAndSolve(config->getSolverFinalConfig(), alpha);
}

std::unique_ptr<RefinementFunctor>
ModelFittingRelativeDensityRatioEstimation::getRefinementFunctor() {
  sgpp::base::AdaptivityConfiguration &refinementConfig = this->config->getRefinementConfig();
  switch (refinementConfig.refinementFunctorType_) {
    case RefinementFunctorType::Surplus: {
      return std::make_unique<SurplusRefinementFunctor>(
          alpha, config->getRefinementConfig().numRefinementPoints_,
          config->getRefinementConfig().refinementThreshold_);
    }
    case RefinementFunctorType::SurplusVolume: {
      return std::make_unique<SurplusVolumeRefinementFunctor>(
          alpha, config->getRefinementConfig().numRefinementPoints_,
          config->getRefinementConfig().refinementThreshold_);
    }
    case RefinementFunctorType::SurplusAbsoluteValue: {
      return std::make_unique<SurplusAbsValueRefinementFunctor>(
          *grid, alpha, config->getRefinementConfig().numRefinementPoints_,
          config->getRefinementConfig().refinementThreshold_);
    }
    case RefinementFunctorType::DataBased: {
      std::string errorMessage =
          "Unsupported refinement functor type DataBased "
          "for relative density ratio estimation!";
      throw application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::ZeroCrossing: {
      std::string errorMessage =
          "Unsupported refinement functor type ZeroCrossing "
          "for relative density ratio estimation!";
      throw application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::MultipleClass: {
      std::string errorMessage =
          "Unsupported refinement functor type MultipleClass "
          "for relative density ratio estimation!";
      throw application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::Classification: {
      std::string errorMessage =
          "Unsupported refinement functor type Classification "
          "for relative density ratio estimation!";
      throw application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::GridPointBased: {
      std::string errorMessage =
          "Unsupported refinement functor type GridPointBased "
          "for relative density ratio estimation!";
      throw application_exception(errorMessage.c_str());
    }
  }
  return nullptr;
}

std::unique_ptr<CoarseningFunctor>
ModelFittingRelativeDensityRatioEstimation::getCoarseningFunctor() {
  sgpp::base::AdaptivityConfiguration &adaptivityConfig = this->config->getRefinementConfig();
  switch (adaptivityConfig.coarseningFunctorType_) {
    case CoarseningFunctorType::Surplus: {
      return std::make_unique<SurplusCoarseningFunctor>(
          alpha, config->getRefinementConfig().numCoarseningPoints_,
          config->getRefinementConfig().coarseningThreshold_);
    }
    case CoarseningFunctorType::SurplusVolume: {
      return std::make_unique<SurplusVolumeCoarseningFunctor>(
          alpha, config->getRefinementConfig().numCoarseningPoints_,
          config->getRefinementConfig().coarseningThreshold_);
    }
    case CoarseningFunctorType::SurplusAbsoluteValue: {
      return std::make_unique<SurplusAbsValueCoarseningFunctor>(
          *grid, alpha, config->getRefinementConfig().numCoarseningPoints_,
          config->getRefinementConfig().coarseningThreshold_);
    }
    case CoarseningFunctorType::Classification: {
      std::string errorMessage =
          "Unsupported coarsening functor type Classification "
          "for relative density ratio estimation!";
      throw application_exception(errorMessage.c_str());
    }
  }
  return nullptr;
}

bool ModelFittingRelativeDensityRatioEstimation::adapt() {
  if (grid != nullptr) {
    if (this->initialGridSize == 0) {
      this->initialGridSize = grid->getSize();
    }

    if (refinementsPerformed < config->getRefinementConfig().numRefinements_) {
      auto oldNoPoints = grid->getSize();
      std::vector<size_t> deletedGridPoints;

      // do coarsening before refinement to prevent refined grid points from being coarsened
      // immediately

      // create coarsening functor
      std::unique_ptr<CoarseningFunctor> coarseningFunc = getCoarseningFunctor();

      if (coarseningFunc) {
        // coarsen grid
        if (config->getRefinementConfig().coarsenInitialPoints_) {
          grid->getGenerator().coarsenNFirstOnly(*coarseningFunc, grid->getSize(),
                                                 &deletedGridPoints, 0);
        } else {
          grid->getGenerator().coarsenNFirstOnly(*coarseningFunc, grid->getSize(),
                                                 &deletedGridPoints, this->initialGridSize);
        }
      } else {
        throw application_exception(
            "ModelFittingRelativeDensityRatioEstimation: No coarsening functor could be "
            "created!");
      }

      // create refinemwnt functor (to get updated/coarsened grid!!!)
      std::unique_ptr<RefinementFunctor> refinementFunc = getRefinementFunctor();

      if (refinementFunc) {
        // refine grid
        std::cout << "Old number points " << oldNoPoints << std::endl;
        GeometryConfiguration geometryConfig = config->getGeometryConfig();
        if (!geometryConfig.stencils_.empty()) {
          GridFactory gridFactory;
          grid->getGenerator().refineInter(*refinementFunc,
                                           gridFactory.getInteractions(geometryConfig));
        } else {
          grid->getGenerator().refine(*refinementFunc);
        }
      } else {
        throw application_exception(
            "ModelFittingRelativeDensityRatioEstimation: No refinement functor could be "
            "created!");
      }

      auto newNoPoints = grid->getSize();
      std::cout << "New number points " << newNoPoints << std::endl;
      if (newNoPoints != oldNoPoints || !deletedGridPoints.empty()) {
        // Tell the SLE manager that the grid changed (for internal data structures)

        // Coarsening, remove idx from alpha
        if (deletedGridPoints.size() > 0) {
          // Restructure alpha
          alpha.remove(deletedGridPoints);
        }

        // oldNoPoint refers to the grid size after coarsening
        auto oldNoPoints2 = alpha.size();

        // Refinement, expand alpha
        if (newNoPoints > oldNoPoints2) {
          alpha.resizeZero(grid->getSize());
        }

        // Solve on new grid
        assembleSystemAndSolve(config->getSolverRefineConfig(), alpha);

        refinementsPerformed++;
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  } else {
    throw application_exception(
        "ModelFittingDensityRatioEstimation: Can't adapt before initial grid is created");
  }
}

void ModelFittingRelativeDensityRatioEstimation::update(Dataset &newDatasetP,
                                                        Dataset &newDatasetQ) {
  if (grid != nullptr) {
    // reassign datasets
    dataset = &newDatasetP;
    extraDataset = &newDatasetQ;
    // create sytem matrix
    assembleSystemAndSolve(config->getSolverFinalConfig(), alpha);
  } else {
    fit(newDatasetP, newDatasetQ);
  }
}

DMSystemMatrixTwoDatasets *ModelFittingRelativeDensityRatioEstimation::buildSystemMatrix(
    Grid &grid, DataMatrix &trainDatasetP, DataMatrix &trainDatasetQ, double lambda,
    OperationMultipleEvalConfiguration &mutipleEvalconfig) const {
  auto systemMatrix = new SystemMatrixRelativeDensityRatioEstimation(
      grid, trainDatasetP, trainDatasetQ, lambda, config->getDensityEstimationConfig().omega_);
  systemMatrix->setImplementation(mutipleEvalconfig);

  return systemMatrix;
}

void ModelFittingRelativeDensityRatioEstimation::reset() {
  grid.reset();
  refinementsPerformed = 0;
}

void ModelFittingRelativeDensityRatioEstimation::assembleSystemAndSolve(
    const SLESolverConfiguration &solverConfig, DataVector &alpha) const {
  auto systemMatrix = std::unique_ptr<DMSystemMatrixTwoDatasets>(buildSystemMatrix(
      *grid, dataset->getData(), extraDataset->getData(), config->getRegularizationConfig().lambda_,
      config->getMultipleEvalConfig()));

  DataVector b(grid->getSize());
  systemMatrix->generateb(b);

  reconfigureSolver(*solver, solverConfig);
  solver->solve(*systemMatrix, alpha, b, true, verboseSolver, DEFAULT_RES_THRESHOLD);
}
}  // namespace datadriven
}  // namespace sgpp
