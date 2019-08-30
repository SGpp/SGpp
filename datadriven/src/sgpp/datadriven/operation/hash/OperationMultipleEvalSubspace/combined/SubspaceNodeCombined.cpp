// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <limits>

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/combined/OperationMultipleEvalSubspaceCombinedParameters.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/combined/SubspaceNodeCombined.hpp>

namespace sgpp {
namespace datadriven {

SubspaceNodeCombined::SubspaceNodeCombined(std::vector<uint32_t>& level, uint32_t flatLevel,
                                           std::vector<uint32_t>& hInverse,
                                           std::vector<uint32_t>& index) {
  size_t dim = level.size();
  this->level = level;
  this->hInverse = hInverse;
  this->indices = index;

  // exactly one gp is added in the loop above
  this->existingGridPointsOnLevel = 1;
  this->flatLevel = flatLevel;
  this->type = NOT_SET;

  this->gridPointsOnLevel = 1;

  for (size_t j = 0; j < dim; j++) {
    uint32_t dimTemp = hInverse[j];
    dimTemp >>= 1;  // skip even indices
    this->gridPointsOnLevel *= dimTemp;
  }

  // initalize other member variables with dummies
  this->jumpTargetIndex = 9999;
  this->arriveDiff = 9999;

  // initialize the lock for this subspace
  omp_init_lock(&this->subspaceLock);
}

SubspaceNodeCombined::SubspaceNodeCombined(size_t dim, uint32_t index) {
  for (size_t i = 0; i < dim; i++) {
    this->level.push_back(1);
    this->hInverse.push_back(2);
  }

  this->gridPointsOnLevel = 0;
  this->existingGridPointsOnLevel = 0;
  this->flatLevel = 0;
  this->type = NOT_SET;

  // initalize other member variables with dummies
  this->jumpTargetIndex = index;
  this->arriveDiff = 9999;
}

void SubspaceNodeCombined::lockSubspace() { omp_set_lock(&this->subspaceLock); }

void SubspaceNodeCombined::unlockSubspace() { omp_unset_lock(&this->subspaceLock); }

// increases number of grid points on the subspace
void SubspaceNodeCombined::addGridPoint(std::vector<uint32_t>& index) {
  size_t dim = index.size();

  for (size_t i = 0; i < dim; i++) {
    this->indices.push_back(index[i]);
  }

  this->existingGridPointsOnLevel += 1;
}

void SubspaceNodeCombined::printLevel() {
  for (size_t i = 0; i < level.size(); i++) {
    if (i > 0) {
      std::cout << ", ";
    }

    std::cout << level[i];
  }

  std::cout << std::endl;
}

// unpack has to be called when the subspace is set up (except for surplus valus)
// this method will decide how to best represent the subspace (list or array type)
// and prepare the subspace for its representation
void SubspaceNodeCombined::unpack() {
  double usageRatio = (double)this->existingGridPointsOnLevel / (double)this->gridPointsOnLevel;

  if (usageRatio < X86COMBINED_LIST_RATIO &&
      this->existingGridPointsOnLevel < X86COMBINED_STREAMING_THRESHOLD) {
    this->type = LIST;
  } else {
    this->type = ARRAY;

    this->subspaceArray.resize(this->gridPointsOnLevel);

    for (size_t i = 0; i < this->gridPointsOnLevel; i++) {
      this->subspaceArray[i] = std::numeric_limits<double>::quiet_NaN();
    }
  }
}

// the first call initializes the array for ARRAY type subspaces
//
void SubspaceNodeCombined::setSurplus(size_t indexFlat, double surplus) {
  if (this->type == ARRAY) {
    this->subspaceArray[indexFlat] = surplus;
  } else if (this->type == LIST) {
    bool found = false;

    for (std::pair<uint32_t, double>& tuple : this->indexFlatSurplusPairs) {
      if (tuple.first == indexFlat) {
        tuple.second = surplus;
        found = true;
        break;
      }
    }

    if (!found) {
      this->indexFlatSurplusPairs.emplace_back(std::make_pair(indexFlat, surplus));
    }
  }
}

// the first call initializes the array for ARRAY type subspaces
double SubspaceNodeCombined::getSurplus(size_t indexFlat) {
  if (this->type == ARRAY) {
    return this->subspaceArray[indexFlat];
  } else if (this->type == LIST) {
    for (std::pair<uint32_t, double> tuple : this->indexFlatSurplusPairs) {
      if (tuple.first == indexFlat) {
        return tuple.second;
      }
    }
  }

  throw;
}

uint32_t SubspaceNodeCombined::compareLexicographically(SubspaceNodeCombined& current,
                                                        SubspaceNodeCombined& last) {
  for (uint32_t i = 0; i < current.level.size(); i++) {
    if (current.level[i] != last.level[i]) {
      return i;
    }
  }

  throw "illegal input";
}

bool SubspaceNodeCombined::subspaceCompare(SubspaceNodeCombined left, SubspaceNodeCombined right) {
  for (size_t i = 0; i < left.level.size(); i++) {
    if (left.level[i] >= right.level[i]) {
      if (left.level[i] > right.level[i]) {
        return 0;
      }
    } else {
      return 1;
    }
  }

  return 1;
}
}
}
