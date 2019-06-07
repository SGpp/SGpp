/**
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * GridFactory.hpp
 *
 *  Created on: May 22, 2018
 *      Author: dominik
 */

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/configuration/GeometryConfiguration.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Factory class to create grids based on a configuration file
 */
class GridFactory {
 public:
  /*
   * Default constructor
   */
  GridFactory() {}

  /**
   * Creates and initializes a grid based on a configuration file
   * @param gridConfig the grid configuration
   * @param interactions the interactions for each dimension
   * @return pointer to the grid object created
   */
  sgpp::base::Grid* createGrid(const sgpp::base::GeneralGridConfiguration& gridConfig,
                               const std::vector<std::vector<size_t>> interactions) const;

  /*
   * method to decide which interactions have to be calculated based on the stencil
   * @param stencil(geometry relation of pixels) e.g. DirectNeighbours
   * @return returns the calculated interaction that have been choosen by the stencil
   */
  std::vector<std::vector<size_t>> getInteractions(
      sgpp::datadriven::StencilType stencilType,
      std::vector<std::vector<int64_t>>& imageDimensions) const;

  /*
   * calculates hierachical parent interactions
   * @param vector of resolution
   * @return all hierachical parent  interactions of all pixels in a vector
   */
  std::vector<std::vector<size_t>> getHierachicalParents(
      std::vector<std::vector<int64_t>>& imageDimensions) const;

  /*
   * calculates direct neighbour interactions
   * @param vector of resolution
   * @return all direct neighbour interactions of all pixels in a vector
   */
  std::vector<std::vector<size_t>> getDirectNeighbours(
      std::vector<std::vector<int64_t>>& imageDimensions) const;

 private:
  void addChildParentInteractionRecursive(
      std::vector<double>* rescale, std::vector<int64_t>* childDim, size_t currentDimension,
      std::vector<int64_t>* childPosition, std::vector<int64_t>* parentPosition,
      std::vector<size_t>* parentMultiplicators, std::vector<size_t>* childMultiplicators,
      size_t parentOffset, size_t childOffset, std::vector<std::vector<size_t>>* res) const;

  
  /*
   * calculates the index of a given position. This method exspect a row wise data layout
   * with each additional dimension beeing stored between the data of the pixels
   * @param vector containing the size of each dimension
   * @param vector containing the position
   * @return the index of the given position in the data
   */
  size_t getDataIndex(size_t numberOfDimensions, std::vector<size_t>* multiplicators,
                      std::vector<int64_t>* position) const;
  size_t getMultiplicator(std::vector<int64_t> imageDimensions, size_t dimension) const;
  void getNextPosition(std::vector<int64_t>* dimension, std::vector<int64_t>* position) const;
  
};
}  // namespace datadriven
}  // namespace sgpp
