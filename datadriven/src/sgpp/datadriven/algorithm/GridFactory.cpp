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
#include <sgpp/base/grid/GridDataBase.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>

#include <string>
#include <vector>

using sgpp::base::algorithm_exception;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridType;

namespace sgpp {
namespace datadriven {

sgpp::base::Grid* GridFactory::createGrid(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const std::vector<std::vector<size_t>> interactions) const {
  Grid* tmpGrid;
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
    throw algorithm_exception("LearnerBase::InitializeGrid: An unsupported grid type was chosen!");
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
    sgpp::datadriven::StencilType stencilType, std::vector<std::vector<int64_t>>& dim) const {
  std::vector<std::vector<size_t>> interactions;

  std::vector<std::vector<int64_t>> res = dim;

  switch (stencilType) {
    case sgpp::datadriven::StencilType::DN:
      interactions = getDirectNeighbours(res);
      break;
    case sgpp::datadriven::StencilType::HP:
      interactions = getHierachicalParents(res);
      break;
    case sgpp::datadriven::StencilType::None:
      interactions = std::vector<std::vector<size_t>>();
      break;
    default:
      std::cout << "Stencil not found";
      std::cout << std::endl;
      break;
  }

  return interactions;
}

std::vector<std::vector<size_t>> sgpp::datadriven::GridFactory::getHierachicalParents(
    std::vector<std::vector<int64_t>>& res) const {
  std::vector<std::vector<size_t>> vec = std::vector<std::vector<size_t>>();
  size_t childOffset = 0;
  std::vector<size_t> childMultiplicators = std::vector<size_t>();
  std::vector<int64_t> childDim;

  for (auto dimension : res) {
    std::vector<size_t> multiplicators = std::vector<size_t>();
    multiplicators.push_back(1);
    for (size_t i = 1; i < dimension.size(); i++) {
      multiplicators.push_back(multiplicators.at(i - 1) * dimension.at(i - 1));
    }

    size_t offset = 0;

    // skip first image
    if (childMultiplicators.size() != 0) {
      // Prepare some values
      offset = childOffset + childMultiplicators.at(childMultiplicators.size() - 1) *
                                 childDim.at(childDim.size() - 1);
      std::vector<double> rescaled = std::vector<double>();
      for (size_t i = 0; i < dimension.size(); i++) {
        rescaled.push_back(static_cast<double>(dimension.at(i)) /
                           static_cast<double>(childDim.at(i)));
      }

      std::vector<int64_t> position = std::vector<int64_t>(dimension.size(), 0);
      size_t dataDimension = 1;
      for (size_t i = 0; i < dimension.size(); i++) {
        dataDimension *= dimension.at(i);
      }

      // Iterate over all cells and add its possible parents
      do {
        std::vector<int64_t> parentPosition = std::vector<int64_t>(dimension.size(), 0);
        addChildParentInteractionRecursive(&rescaled, &childDim, 0, &position, &parentPosition,
                                           &multiplicators, &childMultiplicators, offset,
                                           childOffset, &vec);

        getNextPosition(&childDim, &position);
      } while (position.size() != 0);
    }

    // 1d vector for all dimensions
    for (size_t i = 0;
         i < multiplicators.at(dimension.size() - 1) * dimension.at(dimension.size() - 1); i++) {
      std::vector<size_t> tmp = std::vector<size_t>();
      tmp.push_back(offset + i);
      vec.push_back(tmp);
    }

    childMultiplicators = multiplicators;
    childDim = dimension;
    childOffset = offset;
  }

  // add empty vector
  std::vector<size_t> empty = std::vector<size_t>();
  vec.push_back(empty);

  return vec;
}

void sgpp::datadriven::GridFactory::addChildParentInteractionRecursive(
    std::vector<double>* rescale, std::vector<int64_t>* childDim, size_t currentDimension,
    std::vector<int64_t>* childPosition, std::vector<int64_t>* parentPosition,
    std::vector<size_t>* parentMultiplicators, std::vector<size_t>* childMultiplicators,
    size_t parentOffset, size_t childOffset, std::vector<std::vector<size_t>>* res) const {
  if (currentDimension < childDim->size()) {
    size_t position = static_cast<size_t>(rescale->at(currentDimension) *
                                          static_cast<double>(childPosition->at(currentDimension)));
    double nextPosition = rescale->at(currentDimension) *
                          static_cast<double>(childPosition->at(currentDimension) + 1);
    // iterate over all positions inbetween (neccesarry for increasing dimensions)
    for (size_t i = position; static_cast<double>(i) < nextPosition; i++) {
      parentPosition->at(currentDimension) = i;
      addChildParentInteractionRecursive(rescale, childDim, currentDimension + 1, childPosition,
                                         parentPosition, parentMultiplicators, childMultiplicators,
                                         parentOffset, childOffset, res);
    }
  } else {
    auto tmp = std::vector<size_t>();
    tmp.push_back(getDataIndex(childDim->size(), childMultiplicators, childPosition) + childOffset);
    tmp.push_back(getDataIndex(childDim->size(), parentMultiplicators, parentPosition) +
                  parentOffset);
    res->push_back(tmp);
  }
}

size_t sgpp::datadriven::GridFactory::getDataIndex(size_t numberOfDimensions,
                                                   std::vector<size_t>* multiplicators,
                                                   std::vector<int64_t>* position) const {
  size_t index = 0;
  for (size_t i = 0; i < numberOfDimensions; i++) {
    index += multiplicators->at(i) * position->at(i);
  }
  return index;
}

void sgpp::datadriven::GridFactory::getNextPosition(std::vector<int64_t>* dimension,
                                                    std::vector<int64_t>* position) const {
  for (size_t i = 0; i < position->size(); i++) {
    if (position->at(i) + 1 < dimension->at(i)) {
      position->at(i)++;
      break;
    } else {
      if (i == position->size() - 1)  // last element reached
        position->clear();
      else
        position->at(i) = 0;
    }
  }
}

std::vector<std::vector<size_t>> sgpp::datadriven::GridFactory::getDirectNeighbours(
    std::vector<std::vector<int64_t>>& imageResolutions) const {
  std::vector<std::vector<size_t>> vec = std::vector<std::vector<size_t>>();
  size_t imageOffset = 0;
  for (std::vector<int64_t> res : imageResolutions) {
    std::vector<size_t> multiplicators = std::vector<size_t>();
    multiplicators.push_back(1);
    for (size_t i = 1; i < res.size(); i++) {
      multiplicators.push_back(multiplicators.at(i - 1) * res.at(i - 1));
    }
    std::vector<int64_t> position = std::vector<int64_t>(res.size(), 0);

    do {
      for (size_t i = 0; i < res.size(); i++) {
        if (position.at(i) + 1 < res.at(i)) {
          std::vector<size_t> tmp = std::vector<size_t>();
          tmp.push_back(imageOffset + getDataIndex(res.size(), &multiplicators, &position));
          position.at(i)++;
          tmp.push_back(imageOffset + getDataIndex(res.size(), &multiplicators, &position));
          position.at(i)--;
          vec.push_back(tmp);
        }
      }

      getNextPosition(&res, &position);
    } while (position.size() != 0);

    // 1d vector for all dimensions
    for (size_t i = 0; i < multiplicators.at(res.size() - 1) * res.at(res.size() - 1); i++) {
      std::vector<size_t> tmp = std::vector<size_t>();
      tmp.push_back(imageOffset + i);
      vec.push_back(tmp);
    }

    for (int64_t size : res) {
      imageOffset += size;
    }
  }

  // add empty vector
  std::vector<size_t> empty = std::vector<size_t>();
  vec.push_back(empty);
  return vec;
}
}  // namespace datadriven
}  // namespace sgpp
