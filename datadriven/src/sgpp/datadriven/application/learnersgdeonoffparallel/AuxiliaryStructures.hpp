// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <vector>
#include <list>

namespace sgpp {
namespace datadriven {

/**
 * Compile time option for activating debug messages
 */
#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif


/**
 * Pair to hold the level and the index of a grid point in one dimension
 */
struct LevelIndexPair {
  /**
   * Level of the grid point in the hierarchy.
   */
  sgpp::base::HashGridPoint::level_type level;
  /**
   * Index of the grid point in the index set for that particular level.
   */
  sgpp::base::HashGridPoint::index_type index;
};

/**
 * Vector that holds level index pairs for every dimensions
 */
typedef std::vector<LevelIndexPair> LevelIndexVector;

/**
 * Structure to hold the grid modifications for a refinement cycle for one class
 */
struct RefinementResult {
  /**
   * Holds a list of the positions in hierarchical form for every grid point that was added.
   * These are found at the end of the grid structure after a refinement cycle.
   */
  std::list<LevelIndexVector> addedGridPoints;

  /**
   * Holds the index in the grid interpreted as a vector of a point that was deleted.
   * These should be applied in order.
   */
  std::list<size_t> deletedGridPointsIndices;
};
}  // namespace datadriven
}  // namespace sgpp
