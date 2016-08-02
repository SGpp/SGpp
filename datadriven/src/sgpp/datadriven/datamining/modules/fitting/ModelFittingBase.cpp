/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 *
 *
 * Author: Michael Lettrich
 */

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

namespace sgpp {
namespace datadriven {

using base::DataMatrix;
using base::DataVector;
using base::OperationMatrix;
using base::GridType;
using base::Grid;

using base::GridPoint;
using base::application_exception;
using base::GridGenerator;

ModelFittingBase::ModelFittingBase() : grid(nullptr), alpha(nullptr) {}

ModelFittingBase::~ModelFittingBase() {}

double ModelFittingBase::evaluate(DataVector& sample) {
  // TODO(lettrich): write test

  op_factory::createOperationHierarchisation(*grid)->doHierarchisation(*alpha);
  auto opEval(op_factory::createOperationEval(*grid));
  return opEval->eval(*alpha, sample);
}

std::unique_ptr<DataVector> ModelFittingBase::evaluate(DataMatrix& samples) {
  // TODO(lettrich): write test
  op_factory::createOperationHierarchisation(*grid)->doHierarchisation(*alpha);
  auto opMultEval(op_factory::createOperationMultipleEval(*grid, samples));
  auto result = std::make_unique<DataVector>(samples.getNrows());
  opMultEval->eval(*alpha, *result);
  return result;
}

std::shared_ptr<OperationMatrix> ModelFittingBase::getRegularizationMatrix(
    sgpp::datadriven::RegularizationType regType) {
  std::shared_ptr<OperationMatrix> C;

  if (regType == sgpp::datadriven::RegularizationType::Identity) {
    C = std::shared_ptr<OperationMatrix>(sgpp::op_factory::createOperationIdentity(*grid));
  } else if (regType == sgpp::datadriven::RegularizationType::Laplace) {
    C = std::shared_ptr<OperationMatrix>(sgpp::op_factory::createOperationLaplace(*grid));
  } else {
    throw base::application_exception(
        "ModelFittingBase::getRegularizationMatrix - unknown regularization type");
  }

  return C;
}

void ModelFittingBase::initializeGrid(base::RegularGridConfiguration gridConfig) {
  // load grid
  if (gridConfig.type_ == GridType::Linear) {
    grid = std::shared_ptr<Grid>(Grid::createLinearGrid(gridConfig.dim_));
  } else if (gridConfig.type_ == GridType::LinearL0Boundary) {
    grid = std::shared_ptr<Grid>(Grid::createLinearBoundaryGrid(
        gridConfig.dim_, static_cast<GridPoint::level_type>(gridConfig.boundaryLevel_)));
  } else if (gridConfig.type_ == GridType::LinearBoundary) {
    grid = std::shared_ptr<Grid>(Grid::createLinearBoundaryGrid(gridConfig.dim_));
  } else {
    throw application_exception("ModelFittingBase::createRegularGrid: grid type is not supported");
  }

  GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(gridConfig.level_);
}

std::shared_ptr<sgpp::base::Grid> ModelFittingBase::getGrid() { return grid; }

std::shared_ptr<sgpp::base::DataVector> ModelFittingBase::getSurpluses() { return alpha; }

} /* namespace datadriven */
} /* namespace sgpp */
