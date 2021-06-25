// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <map>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>
#include <cstddef>

namespace sgpp {
namespace datadriven {

/**
 * Hash to make an unordered_set over an vector of size_t's possible.
 */
struct VectorHash {
    size_t operator()(const std::vector<size_t>& v) const {
        std::hash<size_t> hasher;
        size_t seed = 0;
        for (size_t i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

class CombiScheme {
 public:
  /**
   * Empty Constructor
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
   * @return whether we can refine the scheme
   */
  bool isRefinable(std::vector<size_t> levelvec);

  /**
   * update the scheme by component of levelvec
   * @param levelvec vector for the component in question
   * @return whether the scheme was refined
   */
  bool refineComponent(std::vector<size_t> levelvec);

  /**
   * return the combigrid scheme (levelvectors of the components and their coefficients)
   * @return levelvectors of the components and their corresponding coefficients
   */
  std::vector<std::pair<std::vector<size_t>, int>> getCombiScheme();

 private:
  /**
   * Set containing the indices
   */
  std::map<std::vector<size_t>, bool> index_set;

  /**
   * dimension
   */
  size_t dimension = 0;

  /**
   * level
   */
  size_t level = 0;

  /**
   * Initialize the index set
   */
  void initIndexSet();

  /**
   * This method computes recursively all possible level vectors of dimension
   * dim that have an l_1 norm of values.
   * This is used to efficiently compute the standard combination scheme.
   * @param dim dimension
   * @param values values that the l_1 norm of levels vectors of components have to fulfill
   * @return parts to build the index set
   */
  std::unordered_set<std::vector<size_t>, VectorHash> getGrids(size_t dim, size_t values);

  /**
   * update the scheme by component of levelvec
   * @param dim dimension in which to refine
   * @param levelvec vector for the component that should be refined
   */
  void refineScheme(size_t dim, std::vector<size_t> levelvec);
};
} /* namespace datadriven */
} /* namespace sgpp */
