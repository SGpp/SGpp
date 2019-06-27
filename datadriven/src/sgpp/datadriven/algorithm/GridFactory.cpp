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
    sgpp::datadriven::StencilType stencilType, std::vector<std::vector<int64_t>>& dim) const {
  std::vector<std::vector<size_t>> interactions;

  std::vector<std::vector<int64_t>> res = dim;

  switch (stencilType) {
    case sgpp::datadriven::StencilType::DirectNeighbour:
      interactions = getDirectNeighbours(res);
      break;
    case sgpp::datadriven::StencilType::HierarchicalParent:
      interactions = getHierarchicalParents(res, false, true);
      break;
    case sgpp::datadriven::StencilType::RecursiveHierarchicalParent:
      interactions = getHierarchicalParents(res, true, false);
      break;
    case sgpp::datadriven::StencilType::FullyRecursiveHierarchicalParent:
      interactions = getHierarchicalParents(res, false, false);
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

std::vector<std::vector<size_t>> sgpp::datadriven::GridFactory::getHierarchicalParents(
    std::vector<std::vector<int64_t>>& imageDimensions, bool onlyFirstLevel,
    bool onlyFirstParent) const {
  auto vec = std::vector<std::vector<size_t>>();
  auto multiplicatorsPerLevel = std::vector<std::vector<size_t>>();
  for (auto dimension : imageDimensions) {
    auto multiplicators = std::vector<size_t>();
    multiplicators.push_back(1);
    for (size_t i = 1; i < dimension.size(); i++) {
      multiplicators.push_back(multiplicators.at(i - 1) * dimension.at(i - 1));
    }
    multiplicatorsPerLevel.push_back(multiplicators);
  }

  auto offsetsPerLevel = std::vector<size_t>();
  offsetsPerLevel.push_back(0);
  for (size_t i = 1; i < imageDimensions.size(); i++) {
    auto childDimensions = imageDimensions.at(i - 1);
    auto childMultiplicator = multiplicatorsPerLevel.at(i - 1);
    auto numberOfDataColumnsInChild = childMultiplicator.at(childMultiplicator.size() - 1) *
                                      childDimensions.at(childDimensions.size() - 1);
    auto childOffset = offsetsPerLevel.at(i - 1);
    offsetsPerLevel.push_back(childOffset + numberOfDataColumnsInChild);
  }

  for (size_t i = 0; i < imageDimensions.size() - 1; i++) {
    for (size_t j = i + 1; j < imageDimensions.size(); j++) {
      auto scale = std::vector<double>();
      for (size_t k = 0; k < imageDimensions.at(i).size(); k++) {
        scale.push_back(static_cast<double>(imageDimensions.at(j).at(k)) /
                        static_cast<double>(imageDimensions.at(i).at(k)));
      }

      auto position = std::vector<int64_t>(imageDimensions.at(j).size(), 0);

      // Iterate over all cells and add its possible parents
      do {
        std::vector<int64_t> parentPosition = std::vector<int64_t>(position.size(), 0);
        addChildParentInteractionRecursive(&scale, &imageDimensions.at(i), 0, &position,
                                           &parentPosition, &multiplicatorsPerLevel.at(j),
                                           &multiplicatorsPerLevel.at(i), offsetsPerLevel.at(j),
                                           offsetsPerLevel.at(i), &vec);

        getNextPosition(&imageDimensions.at(i), &position);
      } while (position.size() != 0);
      if (onlyFirstParent) {
        break;
      }
    }
    if (onlyFirstLevel) {
      break;
    }
  }

  addOneDimensionalInteractions(&imageDimensions, &vec);
  vec.push_back(std::vector<size_t>());
  return vec;
}

void sgpp::datadriven::GridFactory::addOneDimensionalInteractions(
    std::vector<std::vector<int64_t>>* imageDimensions,
    std::vector<std::vector<size_t>>* vec) const {
  size_t accumulation = 0;
  for (auto resolution : *imageDimensions) {
    size_t factor = 1;
    for (auto dimension : resolution) {
      factor *= dimension;
    }
    accumulation += factor;
  }

  for (size_t i = 0; i < accumulation; i++) {
    vec->push_back(std::vector<size_t>{i});
  }
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

    size_t factor = 1;
    for (int64_t size : res) {
      factor *= size;
    }
    imageOffset += factor;
  }

  addOneDimensionalInteractions(&imageResolutions, &vec);

  // add empty vector
  vec.push_back(std::vector<size_t>());
  return vec;
}
}  // namespace datadriven
}  // namespace sgpp
