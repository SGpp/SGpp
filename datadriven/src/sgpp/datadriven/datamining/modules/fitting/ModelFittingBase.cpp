/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 *
 *
 * Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

namespace sgpp {
namespace datadriven {

using base::DataVector;
using base::GridType;
using base::Grid;
using base::GridPoint;
using base::factory_exception;
using base::GridGenerator;
using base::GeneralGridConfiguration;
using base::application_exception;
using sgpp::solver::SLESolver;
using sgpp::solver::SLESolverType;
using sgpp::solver::ConjugateGradients;
using sgpp::solver::BiCGStab;
using sgpp::solver::SLESolverConfiguration;

ModelFittingBase::ModelFittingBase()
    : verboseSolver{true}, config{nullptr}, dataset{nullptr}, solver{nullptr} {}

const FitterConfiguration &ModelFittingBase::getFitterConfiguration() const { return *config; }

Grid *ModelFittingBase::buildGrid(const GeneralGridConfiguration &gridConfig) const {
  // load grid
  Grid *tmpGrid;
  if (gridConfig.type_ == GridType::Linear) {
    tmpGrid = Grid::createLinearGrid(gridConfig.dim_);
  } else if (gridConfig.type_ == GridType::LinearL0Boundary) {
    tmpGrid = Grid::createLinearBoundaryGrid(
        gridConfig.dim_, static_cast<GridPoint::level_type>(gridConfig.boundaryLevel_));
  } else if (gridConfig.type_ == GridType::LinearBoundary) {
    tmpGrid = Grid::createLinearBoundaryGrid(gridConfig.dim_);
  } else if (gridConfig.type_ == GridType::ModLinear) {
    tmpGrid = Grid::createModLinearGrid(gridConfig.dim_);
 } else {
    throw factory_exception("ModelFittingBase::buildGrid: grid type is not supported");
  }

  GridGenerator &gridGen = tmpGrid->getGenerator();
  if(gridConfig.generalType_ == base::GeneralGridType::RegularSparseGrid){
	  gridGen.regular(gridConfig.level_);
  } else if(gridConfig.generalType_ == base::GeneralGridType::ComponentGrid){
	  gridGen.anisotropicFull(gridConfig.levelVector_);
  }

  return tmpGrid;
}

SLESolver *ModelFittingBase::buildSolver(const SLESolverConfiguration &sleConfig) const {
  if (sleConfig.type_ == SLESolverType::CG) {
    return new ConjugateGradients(sleConfig.maxIterations_, sleConfig.eps_);
  } else if (sleConfig.type_ == SLESolverType::BiCGSTAB) {
    return new BiCGStab(sleConfig.maxIterations_, sleConfig.eps_);
  } else {
    throw factory_exception(
        "ModelFittingBase: An unsupported SLE solver type was "
            "chosen");
  }
}

void ModelFittingBase::reconfigureSolver(SLESolver &solver,
                                         const SLESolverConfiguration &sleConfig) const {
  solver.setMaxIterations(sleConfig.maxIterations_);
  solver.setEpsilon(sleConfig.eps_);
}
} /* namespace datadriven */
} /* namespace sgpp */
