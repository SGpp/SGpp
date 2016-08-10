// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationMakePositive.hpp>

#include <limits>

namespace sgpp {
namespace datadriven {

class OperationLimitFunctionValueRange {
 public:
  OperationLimitFunctionValueRange(
      base::Grid& grid, datadriven::MakePositiveCandidateSearchAlgorithm candiateSearch =
                            MakePositiveCandidateSearchAlgorithm::Intersections,
      datadriven::MakePositiveInterpolationAlgorithm interpolationAlgorithm =
          MakePositiveInterpolationAlgorithm::SetToZero,
      bool generateConsistentGrid = true, bool verbose = false);
  virtual ~OperationLimitFunctionValueRange();

  void doLowerLimitation(base::Grid*& newGrid, base::DataVector& newAlpha, double ylower,
                         bool resetGrid = true);
  void doUpperLimitation(base::Grid*& newGrid, base::DataVector& newAlpha, double yupper,
                         bool resetGrid = true);
  void doLimitation(base::Grid*& newGrid, base::DataVector& newAlpha, double ylower, double yupper);

 private:
  void addConst(base::Grid& grid, base::DataVector& alpha, double c, double y);

  base::Grid& grid;
  datadriven::MakePositiveCandidateSearchAlgorithm candidateSearch;
  datadriven::MakePositiveInterpolationAlgorithm interpolationAlgorithm;
  bool generateConsistentGrid;
  bool verbose;
};

} /* namespace datadriven */
} /* namespace sgpp */
