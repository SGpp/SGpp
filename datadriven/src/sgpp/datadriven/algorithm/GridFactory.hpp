// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/configuration/GeometryConfiguration.hpp>

#include <set>
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
                               const std::set<std::set<size_t>> interactions) const;

  /*
   * method to decide which interactions have to be calculated based on the stencil
   * @param stencil(geometry relation of pixels) e.g. DirectNeighbours
   * @return returns the calculated interaction that have been choosen by the stencil
   */
  std::set<std::set<size_t>> getInteractions(sgpp::datadriven::GeometryConfiguration config) const;

  void getBlockInteractions(std::set<std::set<size_t>>& interactions,
                            sgpp::datadriven::GeometryConfiguration& geometryConfig,
                            sgpp::datadriven::StencilConfiguration& stencilConf) const;

  /*
   * calculates hierachical parent interactions
   * @param vector of resolution
   * @param defines if only for the first level the hierarchical interactions are created
   * @param defines if only for the first parent each the hierarchical interactions are created
   * @return all hierachical interactions of all pixels in a vector
   */
  void getHierarchicalParents(std::set<std::set<size_t>>& interactions,
                              sgpp::datadriven::GeometryConfiguration& geometryConfig,
                              sgpp::datadriven::StencilConfiguration& stencilConf) const;

  /*
   * calculates direct neighbour interactions
   * @param vector of resolution
   * @return all direct neighbour interactions of all pixels in a vector
   */
  void getDirectNeighbours(std::set<std::set<size_t>>& interactions,
                           sgpp::datadriven::GeometryConfiguration& geometryConfig,
                           sgpp::datadriven::StencilConfiguration& stencilConf) const;

 private:
  /*
   * Adds interactions between childs and parent pixels to the set of interactions by
   * calculating and setting the position per axis of the data
   * @param[in] ratio of the axes sizes of the two layers
   * @param[in] number of axes in the data
   * @param[in] current axis the recursive function is on
   * @param[in] position on the child layer
   * @param[in,out] position on the parent layer
   * @param[in] multiplicator of the parent layer
   * @param[in] multiplicator of the child layer
   * @param[in] data offset of the parent layer
   * @param[in] data offset of the child layer
   * @param[in,out] set of interaction terms
   */
  void addChildParentInteraction(std::vector<double>& ratio, size_t numberOfAxes,
                                 size_t currentAxis, std::vector<int64_t>& childPosition,
                                 std::vector<int64_t>& parentPosition,
                                 std::vector<size_t>& parentMultiplicators,
                                 std::vector<size_t>& childMultiplicators, size_t parentOffset,
                                 size_t childOffset,
                                 std::set<std::set<size_t>>& interactions) const;

  /*
   * calculates the index of a given position. This method exspect a row wise data layout
   * with each additional dimension beeing stored between the data of the pixels
   * @param vector containing the size of each dimension
   * @param vector containing the position
   * @return the index of the given position in the data
   */
  size_t getDataIndex(size_t numberOfDimensions, std::vector<size_t>& multiplicators,
                      std::vector<int64_t>& position) const;

  /*
   * Increases a given position like a counter. Therefore if a overflow happens the position will be
   * set to 0 in all dimensions
   * @param vector of resolution
   * @param current position
   */
  void getNextPosition(std::vector<int64_t>& dimension, std::vector<int64_t>& position,
                       size_t colorIndex) const;

  /*
   * Calculates the multiplicators for each level. The multiplicator is the change in the data index
   * for one step further into each dimension
   * @param vector of resolutions
   * @return vector of vectors containing the multiplicator for each dimension
   */
  std::vector<std::vector<size_t>> getMultiplicatorsPerLevel(
      std::vector<std::vector<int64_t>>& imageDimensions) const;

  /*
   * Adds all color interactions in the respective layer to the set of interactions
   * @param[in] vector of resolutions
   * @param[in] index of the color axis
   * @param[in] data offset of the current layer
   * @param[in] vector of the multiplicators of the current layer
   * @param[in,out] set of interaction terms
   */
  void addColorInteractions(std::vector<int64_t>& layerDim, size_t colorIndex, size_t offset,
                            std::vector<size_t>& multiplicators,
                            std::set<std::set<size_t>>& interactions) const;
  /*
   * Calculates the offset for each level. The offset is the starting index of each image level
   * @param vector of resolutions
   * @param multiplicators per level
   * @return vector containing the offset for each image level
   */
  std::vector<size_t> getOffsetPerLevel(
      std::vector<std::vector<int64_t>>& imageDimensions,
      std::vector<std::vector<size_t>>& multiplicatorsPerLevel) const;

  /*
   * Adds all one dimensional interactions to the resulting vector
   * @param vector of resolutions
   * @param the vector in which the interactions are stored
   */
  void addOneDimensionalInteractions(std::vector<int64_t>& layerDimensions, size_t layerOffset,
                                     std::set<std::set<size_t>>& vec) const;
};
}  // namespace datadriven
}  // namespace sgpp
