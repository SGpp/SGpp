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

#include <sgpp/datadriven/algorithm/GridFactory.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/grid/GridDataBase.hpp>

#include <vector>
#include <string>

using sgpp::base::GridType;
using sgpp::base::Grid;
using sgpp::base::GridType;
using sgpp::base::GridGenerator;
using sgpp::base::algorithm_exception;

namespace sgpp {
namespace datadriven {

sgpp::base::Grid *GridFactory::createGrid(const sgpp::base::GeneralGridConfiguration& gridConfig,
    const std::vector<std::vector <size_t>> interactions) const {
  Grid *tmpGrid;
  if (gridConfig.type_ == GridType::Linear) {
    tmpGrid = Grid::createLinearGrid(gridConfig.dim_);
  } else if (gridConfig.type_ == GridType::LinearL0Boundary) {
    tmpGrid = Grid::createLinearBoundaryGrid(
        gridConfig.dim_, static_cast<base::GridPoint::level_type>(gridConfig.boundaryLevel_));
  } else if (gridConfig.type_ == GridType::LinearBoundary) {
    tmpGrid = Grid::createLinearBoundaryGrid(gridConfig.dim_);
  } else if (gridConfig.type_ == GridType::ModLinear) {
    tmpGrid = Grid::createModLinearGrid(gridConfig.dim_);
  } else {
    throw algorithm_exception(
        "LearnerBase::InitializeGrid: An unsupported grid type was chosen!");
  }

  // Generate component grid
  if (gridConfig.generalType_ == sgpp::base::GeneralGridType::ComponentGrid) {
    tmpGrid->getGenerator().anisotropicFull(gridConfig.levelVector_);
    return tmpGrid;
  }

  // Generate regular Grid with LEVELS Levels
  if (interactions.size() == 0) {
    tmpGrid->getGenerator().regular(gridConfig.level_);
  } else {
    std::cout << "Creating geometry aware sparse grid..." << std::endl;
    tmpGrid->getGenerator().regularInter(gridConfig.level_, interactions, 0.0);
    std::cout << "Interactions set!" << std::endl;
  }
  return tmpGrid;
}

std::vector<std::vector<size_t>> sgpp::datadriven::GridFactory::getInteractions(
    sgpp::datadriven::StencilType stencilType, std::vector<int64_t>& dim) const {
  std::vector<std::vector<size_t>> interactions;

  std::vector<int64_t> res = dim;

  if (stencilType == sgpp::datadriven::StencilType::DN) {
    interactions = getDirectNeighbours(res);
  } else {
    std::cout << "Stencil not found";
    std::cout << std::endl;
  }

  return interactions;
}

std::vector<std::vector<size_t>> sgpp::datadriven::GridFactory::getDirectNeighbours(
    std::vector<int64_t>& res) const {
  std::vector<int64_t> geodim1 = res;
  size_t geodimX = res.at(0);
  size_t geodimY = res.at(1);
  std::vector<std::vector<size_t>> vec = std::vector<std::vector<size_t>>();

  for (size_t i = 0; i < geodimY; i++) {
    for (size_t j = 0; j < geodimX-1; j++) {
      std::vector<size_t> xdir = std::vector<size_t>();

      xdir.push_back(i*geodimX+j);

      xdir.push_back(i*geodimX+j+1);

      vec.push_back(xdir);
    }
  }

  for (size_t i = 0; i < geodimX; i++) {
    for (size_t j = 0; j < geodimY-1; j++) {
      std::vector<size_t> ydir = std::vector<size_t>();

      ydir.push_back(i+j*geodimX);

      ydir.push_back(i+(j+1)*geodimX);

      vec.push_back(ydir);
    }
  }

  // 1d vector for all dimensions
  for (size_t i = 0; i < geodimX*geodimY; i++) {
    std::vector<size_t> tmp = std::vector<size_t>();
    tmp.push_back(i);
    vec.push_back(tmp);
  }

  // add empty vector
  std::vector<size_t> empty = std::vector<size_t>();
  vec.push_back(empty);

  return vec;
}
}  // namespace datadriven
}  // namespace sgpp
