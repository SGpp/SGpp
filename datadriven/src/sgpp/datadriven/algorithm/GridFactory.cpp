// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridDataBase.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>

#include <set>
#include <string>
#include <vector>

using sgpp::base::algorithm_exception;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridType;

namespace sgpp {
namespace datadriven {

sgpp::base::Grid* GridFactory::createGrid(const sgpp::base::GeneralGridConfiguration& gridConfig,
                                          const std::set<std::set<size_t>> interactions) const {
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
  std::cout << "Grid Size: " << tmpGrid->getSize() << std::endl;
  return tmpGrid;
}

std::set<std::set<size_t>> sgpp::datadriven::GridFactory::getInteractions(
    sgpp::datadriven::GeometryConfiguration config) const {
  std::set<std::set<size_t>> interactions;
  std::vector<std::vector<int64_t>> res = config.dim_;

  for (StencilConfiguration stencil : config.stencils_) {
    switch (stencil.stencilType_) {
      case sgpp::datadriven::StencilType::DirectNeighbour:
        getDirectNeighbours(interactions, config, stencil);
        break;
      case sgpp::datadriven::StencilType::NextHierarchicalParent:
      case sgpp::datadriven::StencilType::AllHierarchicalParent:
        getHierarchicalParents(interactions, config, stencil);
        break;
      case sgpp::datadriven::StencilType::Block:
        getBlockInteractions(interactions, config, stencil);
        break;
      case sgpp::datadriven::StencilType::None:
        break;
    }
  }
  return interactions;
}

void sgpp::datadriven::GridFactory::getBlockInteractions(
    std::set<std::set<size_t>>& interactions,
    sgpp::datadriven::GeometryConfiguration& geometryConfig,
    sgpp::datadriven::StencilConfiguration& stencilConf) const {
  size_t blockLenght = stencilConf.blockLenght_;
  std::vector<std::vector<size_t>> multiplicatorsPerLevel =
      getMultiplicatorsPerLevel(geometryConfig.dim_);
  std::vector<size_t> offsetsPerLevel =
      getOffsetPerLevel(geometryConfig.dim_, multiplicatorsPerLevel);
  for (size_t i : stencilConf.applyOnLayers_) {
    size_t dimensions = geometryConfig.dim_[i].size();
    size_t offset = offsetsPerLevel[i];
    std::vector<int64_t> position(dimensions, 0);
    while (!position.empty()) {
      for (size_t j = 1; j < static_cast<size_t>(2) << (blockLenght * dimensions); j++) {
        std::set<size_t> tmp;
        for (size_t k = 0; k < blockLenght; k++) {
          for (size_t l = 0; l < dimensions; l++) {
            if (j & static_cast<size_t>(1) << (k * dimensions + l) &&
                static_cast<int64_t>(position[l] + k) < geometryConfig.dim_[i][l]) {
              position[l] += k;
              tmp.insert(getDataIndex(dimensions, multiplicatorsPerLevel[i], position) + offset);
              position[l] -= k;
            }
          }
        }
        interactions.insert(tmp);
      }
      getNextPosition(geometryConfig.dim_[i], position, -1);
    }
    interactions.insert(std::set<size_t>());
  }
}

std::vector<std::vector<size_t>> sgpp::datadriven::GridFactory::getMultiplicatorsPerLevel(
    std::vector<std::vector<int64_t>>& imageDimensions) const {
  std::vector<std::vector<size_t>> multiplicatorsPerLevel;
  for (const std::vector<int64_t>& dimension : imageDimensions) {
    std::vector<size_t> multiplicators;
    multiplicators.push_back(1);
    for (size_t i = 1; i < dimension.size(); i++) {
      multiplicators.push_back(multiplicators[i - 1] * dimension[i - 1]);
    }
    multiplicatorsPerLevel.push_back(multiplicators);
  }
  return multiplicatorsPerLevel;
}

std::vector<size_t> sgpp::datadriven::GridFactory::getOffsetPerLevel(
    std::vector<std::vector<int64_t>>& imageDimensions,
    std::vector<std::vector<size_t>>& multiplicatorsPerLevel) const {
  std::vector<size_t> offsetsPerLevel;
  offsetsPerLevel.push_back(0);
  for (size_t i = 1; i < imageDimensions.size(); i++) {
    std::vector<int64_t> childDimensions = imageDimensions[i - 1];
    std::vector<size_t> childMultiplicator = multiplicatorsPerLevel[i - 1];
    size_t numberOfDataColumnsInChild = childMultiplicator[childMultiplicator.size() - 1] *
                                        childDimensions[childDimensions.size() - 1];
    size_t childOffset = offsetsPerLevel[i - 1];
    offsetsPerLevel.push_back(childOffset + numberOfDataColumnsInChild);
  }
  return offsetsPerLevel;
}

void sgpp::datadriven::GridFactory::getHierarchicalParents(
    std::set<std::set<size_t>>& interactions,
    sgpp::datadriven::GeometryConfiguration& geometryConfig,
    sgpp::datadriven::StencilConfiguration& stencilConf) const {
  std::vector<std::vector<size_t>> multiplicatorsPerLevel =
      getMultiplicatorsPerLevel(geometryConfig.dim_);
  std::vector<size_t> offsetsPerLevel =
      getOffsetPerLevel(geometryConfig.dim_, multiplicatorsPerLevel);

  for (size_t i : stencilConf.applyOnLayers_) {
    for (size_t j = i + 1; j < geometryConfig.dim_.size(); j++) {
      std::vector<double> ratio;
      for (size_t k = 0; k < geometryConfig.dim_[i].size(); k++) {
        ratio.push_back(static_cast<double>(geometryConfig.dim_[j][k]) /
                        static_cast<double>(geometryConfig.dim_[i][k]));
      }

      std::vector<int64_t> position(geometryConfig.dim_[j].size(), 0);

      // Iterate over all positions
      do {
        std::vector<int64_t> parentPosition = std::vector<int64_t>(position.size(), 0);
        addChildParentInteraction(ratio, geometryConfig.dim_[i].size(), 0, position, parentPosition,
                                  multiplicatorsPerLevel[j], multiplicatorsPerLevel[i],
                                  offsetsPerLevel[j], offsetsPerLevel[i], interactions);

        getNextPosition(geometryConfig.dim_[i], position, -1);
      } while (position.size() != 0);

      if (stencilConf.colorIndex_ != -1) {
        addColorInteractions(geometryConfig.dim_[j], stencilConf.colorIndex_, offsetsPerLevel[j],
                             multiplicatorsPerLevel[j], interactions);
      }
      addOneDimensionalInteractions(geometryConfig.dim_[j], offsetsPerLevel[j], interactions);

      if (stencilConf.stencilType_ == StencilType::NextHierarchicalParent) {
        break;
      }
    }
    if (stencilConf.colorIndex_ != -1) {
      addColorInteractions(geometryConfig.dim_[i], stencilConf.colorIndex_, offsetsPerLevel[i],
                           multiplicatorsPerLevel[i], interactions);
    }
    addOneDimensionalInteractions(geometryConfig.dim_[i], offsetsPerLevel[i], interactions);
  }

  // Add empty interaction
  interactions.insert(std::set<size_t>());
}

void sgpp::datadriven::GridFactory::addOneDimensionalInteractions(
    std::vector<int64_t>& imageDimensions, size_t offset, std::set<std::set<size_t>>& vec) const {
  size_t factor = 1;
  for (int64_t dimension : imageDimensions) {
    factor *= dimension;
  }

  for (size_t i = offset; i < factor + offset; i++) {
    vec.insert({i});
  }
}

void sgpp::datadriven::GridFactory::addChildParentInteraction(
    std::vector<double>& ratio, size_t numberOfAxes, size_t currentAxis,
    std::vector<int64_t>& childPosition, std::vector<int64_t>& parentPosition,
    std::vector<size_t>& parentMultiplicators, std::vector<size_t>& childMultiplicators,
    size_t parentOffset, size_t childOffset, std::set<std::set<size_t>>& interactions) const {
  // Check if all axes has been set
  if (currentAxis < numberOfAxes) {
    // Calculate first and last parent on this axis
    size_t begin =
        static_cast<size_t>(ratio[currentAxis] * static_cast<double>(childPosition[currentAxis]));
    double end = ratio[currentAxis] * static_cast<double>(childPosition[currentAxis] + 1);

    // Iterate over all parents' positions
    for (size_t i = begin; static_cast<double>(i) < end; i++) {
      parentPosition[currentAxis] = i;
      addChildParentInteraction(ratio, numberOfAxes, currentAxis + 1, childPosition, parentPosition,
                                parentMultiplicators, childMultiplicators, parentOffset,
                                childOffset, interactions);
    }
  } else {
    std::set<size_t> interaction;
    interaction.insert(getDataIndex(numberOfAxes, childMultiplicators, childPosition) +
                       childOffset);
    interaction.insert(getDataIndex(numberOfAxes, parentMultiplicators, parentPosition) +
                       parentOffset);
    interactions.insert(interaction);
  }
}

size_t sgpp::datadriven::GridFactory::getDataIndex(size_t numberOfDimensions,
                                                   std::vector<size_t>& multiplicators,
                                                   std::vector<int64_t>& position) const {
  size_t index = 0;
  for (size_t i = 0; i < numberOfDimensions; i++) {
    index += multiplicators[i] * position[i];
  }
  return index;
}

void sgpp::datadriven::GridFactory::getNextPosition(std::vector<int64_t>& dimension,
                                                    std::vector<int64_t>& position,
                                                    size_t colorIndex) const {
  for (size_t i = 0; i < position.size(); i++) {
    if (colorIndex == position.size() - 1 && i == colorIndex) {
      position.clear();
    }
    if (i == colorIndex) {
      continue;
    }
    if (position[i] + 1 < dimension[i]) {
      position[i]++;
      break;
    } else {
      if (i == position.size() - 1) {  // last element reached
        position.clear();
      } else {
        position[i] = 0;
      }
    }
  }
}

void sgpp::datadriven::GridFactory::addColorInteractions(
    std::vector<int64_t>& layerDim, size_t colorIndex, size_t offset,
    std::vector<size_t>& multiplicators, std::set<std::set<size_t>>& interactions) const {
  std::vector<int64_t> position = std::vector<int64_t>(layerDim.size(), 0);
  while (position.size() != 0) {
    size_t colorChannels = layerDim[colorIndex];
    for (size_t j = 1; j < (static_cast<size_t>(1) << colorChannels); j++) {
      std::set<size_t> tmp;
      for (size_t k = 0; k < colorChannels; k++) {
        if (j & (1 << k)) {
          position[colorIndex] = k;
          tmp.insert(offset + getDataIndex(layerDim.size(), multiplicators, position));
        }
      }
      interactions.insert(tmp);
    }
    getNextPosition(layerDim, position, colorIndex);
  }
}

void sgpp::datadriven::GridFactory::getDirectNeighbours(
    std::set<std::set<size_t>>& interactions,
    sgpp::datadriven::GeometryConfiguration& geometryConfig,
    sgpp::datadriven::StencilConfiguration& stencilConf) const {
  std::vector<std::vector<size_t>> multiplicatorsPerLevel =
      getMultiplicatorsPerLevel(geometryConfig.dim_);
  std::vector<size_t> offsetsPerLevel =
      getOffsetPerLevel(geometryConfig.dim_, multiplicatorsPerLevel);

  for (size_t i : stencilConf.applyOnLayers_) {
    std::vector<int64_t> position = std::vector<int64_t>(geometryConfig.dim_[i].size(), 0);

    do {
      for (size_t j = 0; j < geometryConfig.dim_[i].size(); j++) {
        if (position[j] + 1 < geometryConfig.dim_[i][j]) {
          std::set<size_t> tmp;
          tmp.insert(offsetsPerLevel[i] + getDataIndex(geometryConfig.dim_[i].size(),
                                                       multiplicatorsPerLevel[i], position));
          position[j]++;
          tmp.insert(offsetsPerLevel[i] + getDataIndex(geometryConfig.dim_[i].size(),
                                                       multiplicatorsPerLevel[i], position));
          position[j]--;
          interactions.insert(tmp);
        }
      }

      // Color index set to -1 to also include neighbours with same color
      getNextPosition(geometryConfig.dim_[i], position, -1);
    } while (position.size() != 0);

    if (stencilConf.colorIndex_ != -1) {
      addColorInteractions(geometryConfig.dim_[i], stencilConf.colorIndex_, offsetsPerLevel[i],
                           multiplicatorsPerLevel[i], interactions);
    }
    addOneDimensionalInteractions(geometryConfig.dim_[i], offsetsPerLevel[i], interactions);
  }

  // add empty vector
  interactions.insert(std::set<size_t>());
}
}  // namespace datadriven
}  // namespace sgpp
