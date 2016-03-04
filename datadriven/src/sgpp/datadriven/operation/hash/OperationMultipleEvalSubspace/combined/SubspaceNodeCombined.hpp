// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <vector>

#include <omp.h>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

class SubspaceNodeCombined {
 public:
  enum SubspaceType {
    NOT_SET, ARRAY, LIST
  };

  std::vector<uint32_t> level;
  std::vector<uint32_t> hInverse;
  uint32_t gridPointsOnLevel;
  uint32_t existingGridPointsOnLevel;
  SubspaceType type;
  //for list representation (and future streaming subspaces)
  std::vector<uint32_t> indices;
  std::vector<std::pair<uint32_t, double> > indexFlatSurplusPairs;
  std::vector<double> subspaceArray;
  omp_lock_t subspaceLock;

  uint32_t jumpTargetIndex;
  uint32_t flatLevel;

  // every node that reaches this subspace has to calculate this diff
  uint32_t arriveDiff;

  SubspaceNodeCombined(std::vector<uint32_t>& level, uint32_t flatLevel,
                       std::vector<uint32_t>& hInverse,
                       std::vector<uint32_t>& index);

  SubspaceNodeCombined(size_t dim, uint32_t index);

  void lockSubspace();

  void unlockSubspace();

  //increases number of grid points on the subspace
  void addGridPoint(std::vector<uint32_t>& index);

  void printLevel();

  // unpack has to be called when the subspace is set up (except for surplus valus)
  // this method will decide how to best represent the subspace (list or array type)
  // and prepare the subspace for its representation
  void unpack();

  // the first call initializes the array for ARRAY type subspaces
  //
  void setSurplus(size_t indexFlat, double surplus);

  // the first call initializes the array for ARRAY type subspaces
  double getSurplus(size_t indexFlat);

  static uint32_t compareLexicographically(SubspaceNodeCombined& current,
      SubspaceNodeCombined& last);

  static bool subspaceCompare(SubspaceNodeCombined left,
                              SubspaceNodeCombined right);

};

}
}
