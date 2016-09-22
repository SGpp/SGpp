// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/tools/NearestNeighbors.hpp>
#include <algorithm>
#include <set>
#include <vector>

namespace sgpp {
namespace datadriven {

NearestNeighbors::NearestNeighbors(size_t rows, size_t cols)
  : rows(rows), cols(cols) {}

NearestNeighbors::point_t NearestNeighbors::idxToPoint(size_t idx) const{
  // Maybe the other way round? :O
  const auto f = idx / rows;
  const auto s = idx % rows;
  return {f, s};
}

size_t NearestNeighbors::pointToIdx(point_t p) const{
    return p.first*rows+p.second;
}

double NearestNeighbors::l2Distance(NearestNeighbors::point_t a,
                                    NearestNeighbors::point_t b) const {
  const auto d1 = a.first - b.first;
  const auto d2 = a.second - b.second;
  return std::sqrt((d1 * d1) + (d2 * d2));
}

std::vector<std::vector<size_t>> NearestNeighbors::getAllInteractions(
    size_t level, double threshold) const {
  const auto neighbors = getAllNeighbors(threshold);
  auto allInteractions = std::vector<std::vector<size_t>>();
  for (size_t i = 1; i < level; ++i) {
    const auto inter = getInteractions(neighbors, i);
    allInteractions.insert(allInteractions.end(), inter.begin(), inter.end());
  }
  // Remove duplicates. This is faster than std::unique.
  auto s = std::set<std::vector<size_t>>(allInteractions.begin(), allInteractions.end());
  allInteractions.assign(s.begin(), s.end());
  return allInteractions;
}

std::vector<size_t> NearestNeighbors::getNeighbors(size_t idx, double threshold) const {
  auto neighbors = std::vector<size_t>();
  // Make sure idx is always first!
  // We need this property later to filter out unwanted interactions.
  neighbors.push_back(idx);
  auto center = idxToPoint(idx);
  for (size_t curIdx = 0; curIdx < rows * cols; ++curIdx) {
    const auto p = idxToPoint(curIdx);
    if (curIdx != idx && l2Distance(center, p) <= threshold) {
      neighbors.push_back(curIdx);
    }
  }
  return neighbors;
}

std::vector<std::vector<size_t>> NearestNeighbors::getAllNeighbors(double threshold) const{
  auto allN = std::vector<std::vector<size_t>>();
  for (size_t i = 0; i < rows * cols; ++i) {
    allN.emplace_back(getNeighbors(i, threshold));
  }
  return allN;
}

std::vector<std::vector<size_t>> NearestNeighbors::getInteractions(
    std::vector<std::vector<size_t>> neighbors, size_t order) const {
  // Follows http://stackoverflow.com/a/9430993/6654913.
  // This function returns all order-long combinations of neighbors.
  // Only combinations that result in points with a valid distance are preserved.
  auto combinations = std::vector<std::vector<size_t>>();
  for (const auto& n : neighbors) {
    auto indices = std::vector<bool>(n.size());
    std::fill(indices.begin(), indices.begin() + order, true);
    const auto first = n[0];  // Only consider interactions where n[0] is contained.
    do {
      auto firstFound = false;
      auto curCombinations = std::vector<size_t>();
      for (size_t i = 0; i < n.size(); ++i) {
        if (indices[i]) {
          curCombinations.push_back(n[i]);
          if (n[i] == first) {
            firstFound = true;
          }
        }
      }
      if (!firstFound) continue;
      // It's easier to remove dumplicates if the vector is sorted.
      std::sort(curCombinations.begin(), curCombinations.end());
      combinations.push_back(curCombinations);
    } while (std::prev_permutation(indices.begin(), indices.end()));
  }
  return combinations;
}
}  // namespace datadriven
}  // namespace sgpp
