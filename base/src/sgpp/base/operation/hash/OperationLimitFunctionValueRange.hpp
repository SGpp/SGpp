// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/operation/hash/OperationMakePositive.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <limits>

namespace sgpp {
namespace base {

class OperationLimitFunctionValueRange {
 public:
  OperationLimitFunctionValueRange(base::Grid& grid,
                                   base::MakePositiveCandidateSearchAlgorithm candiateSearch =
                                       MakePositiveCandidateSearchAlgorithm::Intersections,
                                   base::MakePositiveInterpolationAlgorithm interpolationAlgorithm =
                                       MakePositiveInterpolationAlgorithm::SetToZero,
                                   bool verbose = false);
  virtual ~OperationLimitFunctionValueRange();

  void doLowerLimitation(base::Grid*& newGrid, base::DataVector& newAlpha, double ylower,
                         bool resetGrid = true);
  void doUpperLimitation(base::Grid*& newGrid, base::DataVector& newAlpha, double yupper,
                         bool resetGrid = true);
  void doLimitation(base::Grid*& newGrid, base::DataVector& newAlpha, double ylower, double yupper);

 private:
  void copyGrid(base::Grid& grid, base::Grid*& newGrid);

  void addConst(base::Grid& grid, base::DataVector& alpha, double c, double y);

  base::Grid& grid;
  base::MakePositiveCandidateSearchAlgorithm candidateSearch;
  base::MakePositiveInterpolationAlgorithm interpolationAlgorithm;
  bool verbose;
};

} /* namespace base */
} /* namespace sgpp */
