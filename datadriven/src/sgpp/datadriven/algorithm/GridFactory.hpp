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

#include <vector>
#include <string>

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
  sgpp::base::Grid *createGrid(const sgpp::base::GeneralGridConfiguration& gridConfig,
    const std::vector<std::vector <size_t>> interactions) const;

  /*
   * method to decide which interactions have to be calculated based on the stencil
   * @param stencil(geometry relation of pixels) e.g. DirectNeighbours
   * @return returns the calculated interaction that have been choosen by the stencil
   */
  std::vector<std::vector<size_t>> getInteractions(sgpp::datadriven::StencilType stencilType,
    std::vector<int64_t>& dim) const;

  /*
   * calculates hierachical parent interactions
   * @param vector of resolution
   * @return all hierachical parent  interactions of all pixels in a vector
   */
  std::vector<std::vector<size_t>> getHierachicalParents(std::vector<int64_t>& res) const;

  /*
   * calculates direct neighbour interactions
   * @param vector of resolution
   * @return all direct neighbour interactions of all pixels in a vector
   */
  std::vector<std::vector<size_t>> getDirectNeighbours(std::vector<int64_t>& res) const;

 private:
};
}  // namespace datadriven
}  // namespace sgpp
