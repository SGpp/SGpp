// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NEARESTNEIGHBORS_HPP
#define NEARESTNEIGHBORS_HPP

#include <algorithm>
#include <cmath>
#include <set>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * @brief The NearestNeighbors class
 * @details Generates interactions for images, that only consider
 * interactions between pixels that are close to each other, measured
 * by an \f$ l_2 \f$-distance.
 * The pixel-neighborhood for a threshold of \f$ \sqrt{2} \f$  is
 * \image html nearest_neighbors_l2.svg
 */
class NearestNeighbors {
 public:
  typedef std::pair<size_t, size_t> point_t;

   /**
   * @brief NearestNeighbors
   * @param rows is the number of rows of each image
   * @param cols is the number of columns of each image
   */
  NearestNeighbors(size_t rows, size_t cols);

  /**
   * @brief calculates all interaction terms with l2Distance(center, other) \f$\leq\f$ threshold
   * @param level is the maximum level of interactions to include
   * @param threshold is the maximum distance
   * @return all interactions
   */
  std::vector<std::vector<size_t>> getAllInteractions(size_t level,
                                                      double threshold) const;

 private:
  const size_t rows;
  const size_t cols;
  std::vector<size_t> getNeighbors(size_t idx, double threshold) const;

  std::vector<std::vector<size_t>> getAllNeighbors(double threshold) const;

  std::vector<std::vector<size_t>> getInteractions(
      std::vector<std::vector<size_t>> neighbors, size_t order) const;
  size_t pointToIdx(point_t p) const;
  point_t idxToPoint(size_t idx) const;
  /**
   * @brief Euclidean distance of two pixels
   * @param a first point
   * @param b second point
   * @return the Euclidean distance of a and b
   */
  double l2Distance(point_t a, point_t b) const;
};

}  // namespace datadriven
}  // namespace sgpp
#endif  // NEARESTNEIGHBORS_HPP
