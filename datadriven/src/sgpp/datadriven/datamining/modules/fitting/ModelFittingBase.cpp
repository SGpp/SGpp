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
using base::RegularGridConfiguration;
using base::application_exception;
using sgpp::solver::SLESolver;
using sgpp::solver::SLESolverType;
using sgpp::solver::ConjugateGradients;
using sgpp::solver::BiCGStab;
using sgpp::solver::SLESolverConfiguration;

ModelFittingBase::ModelFittingBase()
    : config{nullptr}, grid{nullptr}, alpha{}, dataset{nullptr}, solver{nullptr} {}

const Grid& ModelFittingBase::getGrid() const {
  if (grid != nullptr) {
    return *grid;
  } else {
    throw application_exception("No grid was fitted yet");
  }
}

const DataVector& ModelFittingBase::getSurpluses() const { return alpha; }

const FitterConfiguration& ModelFittingBase::getFitterConfiguration() const { return *config; }

Grid* ModelFittingBase::buildGrid(const RegularGridConfiguration& gridConfig) const {
  // load grid
  Grid* tmpGrid;
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
    throw factory_exception("ModelFittingBase::createRegularGrid: grid type is not supported");
  }

  GridGenerator& gridGen = tmpGrid->getGenerator();
  gridGen.regular(gridConfig.level_);
  return tmpGrid;
}

SLESolver* ModelFittingBase::buildSolver(const SLESolverConfiguration& sleConfig) const {
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

void ModelFittingBase::reconfigureSolver(SLESolver& solver,
                                         const SLESolverConfiguration& sleConfig) const {
  solver.setMaxIterations(sleConfig.maxIterations_);
  solver.setEpsilon(sleConfig.eps_);
}

} /* namespace datadriven */
} /* namespace sgpp */
