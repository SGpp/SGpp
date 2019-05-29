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
  
  switch (stencilType)
  {
  case sgpp::datadriven::StencilType::DN:
    interactions = getDirectNeighbours(res);
    break;
  case sgpp::datadriven::StencilType::HP:
    interactions = getHierachicalParents(res);
  default:
    std::cout << "Stencil not found";
    std::cout << std::endl;
    break;
  }

  return interactions;
}

std::vector<std::vector<size_t>> sgpp::datadriven::GridFactory::getHierachicalParents(
  std::vector<int64_t>& res) const {
  std::vector<std::vector<size_t>> vec = std::vector<std::vector<size_t>>();
  size_t offset = 0;
  size_t dimX = res.at(0);
  size_t dimY = res.at(1);

  for (size_t i = 2; i < res.size(); i += 2) {
    size_t parentDimX = res.at(i);
    size_t parentDimY = res.at(i + 1);
    size_t parentOffset = dimX * dimY + offset;
    double cellsPerParentX = dimX / static_cast<double>(parentDimX);
    double cellsPerParentY = dimY / static_cast<double>(parentDimY);

    size_t parentY = 0;
    size_t toY = 0;
    double copyCellY = cellsPerParentY;
    for (size_t y = 0; y < dimY; y = toY) {
      toY = ceil(copyCellY) - 1;
      if (toY > dimY) {
        toY = dimY;
      }

      for (size_t j = y; j <= toY; j++) {
        size_t parentX = 0;
        size_t toX = 0;
        double copyCellX = cellsPerParentX;
        for (size_t x = 0; x < dimX; x = toX) {
          toX = ceil(copyCellX) - 1;
          if (toX > dimX) {
            toX = dimX;
          }

          for (size_t i = x; i <= toX; i++) {
            std::vector<size_t> tmp = std::vector<size_t>();
            tmp.push_back(offset + i + dimX * j);
            tmp.push_back(parentOffset + parentX + parentDimX * parentY);
            vec.push_back(tmp);
          }

          if (toX == copyCellX - 1) toX++;
          parentX++;
          copyCellX += cellsPerParentX;
        }
      }
      if (toY == copyCellY - 1) toY++;
      parentY++;
      copyCellY += cellsPerParentY;
    }

    dimX = parentDimX;
    dimY = parentDimY;
    offset = parentOffset;
  }
  return vec;
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
