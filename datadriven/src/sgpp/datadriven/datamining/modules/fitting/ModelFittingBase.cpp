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
    : verboseSolver{true}, config{nullptr}, dataset{nullptr}, solver{nullptr} {}

const FitterConfiguration &ModelFittingBase::getFitterConfiguration() const { return *config; }

Grid *ModelFittingBase::buildGrid(const RegularGridConfiguration &gridConfig) const {
  GridFactory gridFactory;

  // pass interactions with size 0
  std::vector<std::vector <size_t>> interactions = std::vector<std::vector<size_t>>();
  return gridFactory.createGrid(gridConfig, interactions);
}

Grid *ModelFittingBase::buildGrid(const RegularGridConfiguration &gridConfig, const GeometryConfiguration &geometryConfig) const{
  GridFactory gridFactory;
  std::string tmpString = geometryConfig.stencil;
  std::vector<int64_t> dim = geometryConfig.dim;
  return gridFactory.createGrid(gridConfig, gridFactory.getInteractions(tmpString, dim));
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
