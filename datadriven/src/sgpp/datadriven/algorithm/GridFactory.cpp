/**
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * GridFactory.cpp
 *
 *  Created on: May 22, 2018
 *      Author: dominik
 */

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>

#include <vector>

using sgpp::base::GridType;
using sgpp::base::Grid;
using sgpp::base::GridType;
using sgpp::base::GridGenerator;
using sgpp::base::algorithm_exception;

namespace sgpp {
namespace datadriven {

sgpp::base::Grid *GridFactory::createGrid(sgpp::base::GeneralGridConfiguration &gridConfig,
                                          std::vector<std::vector<size_t>> interactions) const {
  Grid *grid;

  if (gridConfig.type_ == GridType::ModLinear) {
    grid = Grid::createModLinearGrid(gridConfig.dim_);
  } else if (gridConfig.type_ == GridType::Linear) {
    grid = Grid::createLinearGrid(gridConfig.dim_);
  } else {
    throw algorithm_exception("LearnerBase::InitializeGrid: An unsupported grid type was chosen!");
  }

  // Generate regular Grid with LEVELS Levels
  /*
  if (interactions.size() == 0) {
    grid->getGenerator().regular(gridConfig.level_);
  } else {
    if (gridConfig.generalType_ != sgpp::base::GeneralGridType::GeometryAwareSparseGrid) {
      // No geometry aware sparse grid but passed interactions nonetheless
      std::cout << "Passed grid configuration for non geometry aware sparse grids and"
                << "interactions vector nonetheless. Geometry aware sparse grid will be created!"
                << std::endl;
    }
    grid->getGenerator().regularInter(gridConfig.level_, interactions, 0.0);
  }
  */

  switch (gridConfig.generalType_) {
    case (sgpp::base::GeneralGridType::RegularSparseGrid):
      if (interactions.size() == 0) {
        grid->getGenerator().regular(gridConfig.level_);
      } else {
        // No geometry aware sparse grid but passed interactions nonetheless
        std::cout << "Passed grid configuration for non geometry aware sparse grids and"
                  << "interactions vector nonetheless. Geometry aware sparse grid will be created!"
                  << std::endl;
        grid->getGenerator().regularInter(gridConfig.level_, interactions, 0.0);
      }
      break;
    case (sgpp::base::GeneralGridType::GeometryAwareSparseGrid):
      if (interactions.size() == 0) {
        grid->getGenerator().regular(gridConfig.level_);
      } else {
        grid->getGenerator().regularInter(gridConfig.level_, interactions, 0.0);
      }
      break;
    case (sgpp::base::GeneralGridType::ComponentGrid):
      grid->getGenerator().anisotropicFull(gridConfig.levelVector_);
      break;
    case (sgpp::base::GeneralGridType::RefinedCoarsenedSparseGrid):
      throw algorithm_exception(
          "GridFactory::createGrid: An unsupported general Grid Type was chosen!");
      break;
    default:
      throw algorithm_exception(
          "GridFactory::createGrid: An unsupported general Grid Type was chosen!");
  }
  return grid;
}
}  // namespace datadriven
}  // namespace sgpp
