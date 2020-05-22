// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <set>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

using base::application_exception;
using base::DataVector;
using base::factory_exception;
using base::GeneralGridConfiguration;
using base::Grid;
using base::GridGenerator;
using base::GridPoint;
using base::GridType;
using sgpp::solver::BiCGStab;
using sgpp::solver::ConjugateGradients;
using sgpp::solver::SLESolver;
using sgpp::solver::SLESolverConfiguration;
using sgpp::solver::SLESolverType;

ModelFittingBase::ModelFittingBase()
    : verboseSolver{true},
      interactions{nullptr},
      config{nullptr},
      dataset{nullptr},
      extraDataset{nullptr},
      solver{nullptr} {}

const FitterConfiguration &ModelFittingBase::getFitterConfiguration() const { return *config; }

FitterConfiguration &ModelFittingBase::getFitterConfiguration() { return *config; }

Grid *ModelFittingBase::buildGrid(const sgpp::base::GeneralGridConfiguration &gridConfig) const {
  GridFactory gridFactory;

  // pass interactions with size 0
  std::set<std::set<size_t>> interactions;
  return gridFactory.createGrid(gridConfig, interactions);
}

Grid *ModelFittingBase::buildGrid(const sgpp::base::GeneralGridConfiguration &gridConfig,
                                  const GeometryConfiguration &geometryConfig) const {
  GridFactory gridFactory;

  // a regular sparse grid is created, if no geometryConfig is defined,
  if (geometryConfig.stencils_.empty()) {
    // interaction with size 0
    std::set<std::set<size_t>> interactions;
    return gridFactory.createGrid(gridConfig, interactions);
  }

  return gridFactory.createGrid(gridConfig, gridFactory.getInteractions(geometryConfig));
}

std::set<std::set<size_t>> ModelFittingBase::getInteractions(
    const GeometryConfiguration &geometryConfig) {
  if (!interactions) {
    if (geometryConfig.stencils_.empty()) {
      // interaction with size 0
      interactions = std::make_unique<std::set<std::set<size_t>>>();
    } else {
      GridFactory gridFactory;
      interactions =
          std::make_unique<std::set<std::set<size_t>>>(gridFactory.getInteractions(geometryConfig));
    }
  }
  return *interactions;
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

Dataset *ModelFittingBase::getDataset() { return dataset; }
} /* namespace datadriven */
} /* namespace sgpp */
