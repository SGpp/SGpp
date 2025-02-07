// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SUBSPACENODESIMPLE_HPP
#define SUBSPACENODESIMPLE_HPP


#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {
class SubspaceNodeSimple {
 public:
  std::vector<size_t> level;
  std::vector<size_t> hInverse;
  size_t actualGridPointsOnLevel;
  std::vector<size_t> indices;  //for stream computations
  size_t gridPointsOnLevel;

  SubspaceNodeSimple(std::vector<size_t>& level, std::vector<size_t>& hInverse,
                     std::vector<size_t>& index) {
    size_t dim = level.size();

    this->level = level;
    this->hInverse = hInverse;
    this->indices = index;

    this->actualGridPointsOnLevel = 1;
    gridPointsOnLevel = 1;

    for (size_t j = 0; j < dim; j++) {
      int dimTemp = static_cast<int>(this->hInverse[j]);
      dimTemp >>= 1;  //skip even indices
      gridPointsOnLevel *= dimTemp;
    }
  }

  // increases number of grid points on the subspace
  void addGridPoint(std::vector<size_t>& index) {
    size_t dim = index.size();

    for (size_t i = 0; i < dim; i++) {
      this->indices.push_back(index[i]);
    }

    this->actualGridPointsOnLevel += 1;
  }
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* SUBSPACENODESIMPLE_HPP */
