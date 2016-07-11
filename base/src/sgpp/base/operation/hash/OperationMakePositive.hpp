// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationMakePositiveCandidateSetAlgorithm.hpp>
#include <sgpp/base/operation/hash/OperationMakePositiveInterpolationAlgorithm.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>
#include <map>

namespace sgpp {
namespace base {

enum class MakePositiveCandidateSearchAlgorithm { FullGrid, Intersections };
enum class MakePositiveInterpolationAlgorithm { SetToZero };

class OperationMakePositive {
 public:
  typedef std::map<size_t, base::HashGridPoint> gridPointCandidatesMap;

  explicit OperationMakePositive(base::Grid& grid,
                                 MakePositiveCandidateSearchAlgorithm candiateSearchAlgorithm =
                                     MakePositiveCandidateSearchAlgorithm::Intersections,
                                 MakePositiveInterpolationAlgorithm interpolationAlgorithm =
                                     MakePositiveInterpolationAlgorithm::SetToZero);

  virtual ~OperationMakePositive();

  /**
   *
   * @param alpha
   * @param newGrod
   * @param newAlpha
   */
  void makePositive(base::Grid*& newGrid, base::DataVector& newAlpha);

 private:
  void copyGrid(base::Grid& grid, base::Grid*& newGrid);

  /**
   *
   * @param alpha
   */
  void makeCurrentNodalValuesPositive(base::DataVector& alpha);

  /**
   *
   * @param alpha
   */
  void forceNewNodalValuesToBePositive(base::Grid& grid, base::DataVector& alpha,
                                       std::vector<size_t>& newGridPoints);

  /**
   *
   * @param candidates
   * @param sortedCandidates
   */
  void sortNonExistingCandidatesByLevelSum(
      std::vector<base::HashGridPoint*>& candidates,
      std::map<size_t, std::vector<base::HashGridPoint*> >& sortedCandidates);

  /**
   *
   * @param alpha
   * @param candidates
   */
  void addFullGridPoints(base::Grid& grid, base::DataVector& alpha,
                         std::vector<base::HashGridPoint*>& candidates,
                         std::vector<size_t>& addedGridPoints);

  /// grid
  base::Grid& grid;

  /// last checked level sum
  size_t lastMinimumCandidateLevelSum;

  /// candidate search algorithm
  std::shared_ptr<base::OperationMakePositiveCandidateSetAlgorithm> candidateSearch;
  std::shared_ptr<base::OperationMakePositiveInterpolationAlgorithm> interpolationMethod;
};

} /* namespace base */
} /* namespace sgpp */
