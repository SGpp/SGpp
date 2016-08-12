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
using base::RegularGridConfiguration;

ModelFittingBase::ModelFittingBase() : grid(nullptr), alpha(nullptr) {}

ModelFittingBase::~ModelFittingBase() {}

double ModelFittingBase::evaluate(DataVector& sample) {
  auto opEval = std::unique_ptr<base::OperationEval>(op_factory::createOperationEval(*grid));
  return opEval->eval(*alpha, sample);
}

void ModelFittingBase::evaluate(DataMatrix& samples, DataVector& results) const {
  auto opMultEval = std::unique_ptr<base::OperationMultipleEval>(
      op_factory::createOperationMultipleEval(*grid, samples));
  opMultEval->eval(*alpha, results);
}

OperationMatrix* ModelFittingBase::getRegularizationMatrix(datadriven::RegularizationType regType) {
  if (regType == datadriven::RegularizationType::Identity) {
    return op_factory::createOperationIdentity(*grid);
  } else if (regType == datadriven::RegularizationType::Laplace) {
    return op_factory::createOperationLaplace(*grid);
  } else {
    throw application_exception(
        "ModelFittingBase::getRegularizationMatrix - unknown regularization type");
  }
}

void ModelFittingBase::initializeGrid(RegularGridConfiguration gridConfig) {
  // load grid
  if (gridConfig.type_ == GridType::Linear) {
    grid = std::shared_ptr<Grid>(Grid::createLinearGrid(gridConfig.dim_));
  } else if (gridConfig.type_ == GridType::LinearL0Boundary) {
    grid = std::shared_ptr<Grid>(Grid::createLinearBoundaryGrid(
        gridConfig.dim_, static_cast<GridPoint::level_type>(gridConfig.boundaryLevel_)));
  } else if (gridConfig.type_ == GridType::LinearBoundary) {
    grid = std::shared_ptr<Grid>(Grid::createLinearBoundaryGrid(gridConfig.dim_));
  } else if (gridConfig.type_ == GridType::ModLinear) {
    grid = std::shared_ptr<Grid>(Grid::createModLinearGrid(gridConfig.dim_));
  } else {
    throw application_exception("ModelFittingBase::createRegularGrid: grid type is not supported");
  }

  GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(gridConfig.level_);
}

const Grid& ModelFittingBase::getGrid() const { return *grid; }

const DataVector& ModelFittingBase::getSurpluses() const { return *alpha; }

} /* namespace datadriven */
} /* namespace sgpp */
