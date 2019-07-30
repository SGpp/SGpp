/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CombiScheme.hpp
 *
 * Created on: Jul 25, 2019
 *     Author: Kilian Röhner
 */

#pragma once

#include <map>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>


namespace sgpp {
namespace datadriven {

class CombiScheme {
 public:
  /**
   * Empty Constructor
   *
   */
  CombiScheme() {}

  /**
   * Initialize the combigrid scheme
   * @param dim dimension
   * @param level level
   */
  void initialize(size_t dim, size_t level);

  /**
   * check if the component of levelvec is refinable
   * @param levelvec vector for the component in question
   * @return wheter we can refine the scheme
   */
  bool isRefinable(std::vector<size_t> levelvec);

  /**
   * update the scheme by component of levelvec
   * @param levelvec vector for the component in question
   * @return wheter the scheme was refined
   */
  bool refineComponent(std::vector<size_t> levelvec);

  /**
   * initialize the active index set
   * @return levelvectors of the components and their corresponding coefficients
   */
  std::vector<std::pair<std::vector<size_t>, int>> getCombiScheme();

 private:
  /**
   * Set containing the active indices
   */
  std::set<std::vector<size_t>> active_index_set;

  /**
   * Set containing the old indices
   */
  std::set<std::vector<size_t>> old_index_set;

  /**
   * dimension
   */
  size_t dimension = 0;

  /**
   * level
   */
  size_t level = 0;

  /**
   * initialize the active index set
   * @return the active index set
   */
  std::set<std::vector<size_t>> init_active_index_set();

  /**
   * initialize the old index set
   * @return the active index set
   */
  std::set<std::vector<size_t>> init_old_index_set();

  /**
   * initialize the active index set
   * @param dim dimension
   * @param values level plane
   * @return the active index set
   */
  std::set<std::vector<size_t>> getGrids(size_t dim, size_t values);

  /**
   * update the scheme by component of levelvec
   * @param dim dimension in which to refine
   * @param levelvec vector for the component in question
   */
  void refine_scheme(size_t dim, std::vector<size_t> levelvec);

  /**
   * initialize the active index set
   * @param grid_dict dictionary containing the component level as key and the coefficients as values
   * @param index_set index set containing the component levels
   */
  void get_coefficients_to_index_set(std::map<std::vector<size_t>, int> &grid_dict, std::set<std::vector<size_t>> &index_set);

};
} /* namespace datadriven */
} /* namespace sgpp */
