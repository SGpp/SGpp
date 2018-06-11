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

#ifndef DATADRIVEN_SRC_SGPP_DATADRIVEN_ALGORITHM_GRIDFACTORY_HPP_
#define DATADRIVEN_SRC_SGPP_DATADRIVEN_ALGORITHM_GRIDFACTORY_HPP_

#include <sgpp/base/grid/Grid.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {


/**
 * Factory class to create grids based on a configuration file
 */
class GridFactory {
 public:
  /**
   * Creates and initializes a grid based on a configuration file
   * @param gridConfig the grid configuration
   * @param interactions the interactions for each dimension
   * @return pointer to the grid object created
   */
  sgpp::base::Grid *createGrid(sgpp::base::GeneralGridConfiguration& gridConfig,
      std::vector<std::vector <size_t>> interactions) const;
};
}  // namespace datadriven
}  // namespace sgpp

#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_ALGORITHM_GRIDFACTORY_HPP_ */
