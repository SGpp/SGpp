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
      datadriven::MakePositiveCandidateSearchAlgorithm candiateSearch =
          MakePositiveCandidateSearchAlgorithm::Intersections,
      datadriven::MakePositiveInterpolationAlgorithm interpolationAlgorithm =
          MakePositiveInterpolationAlgorithm::SetToZero,
      bool verbose = false);

  virtual ~OperationLimitFunctionValueRange();

  void doLowerLimitation(base::Grid& grid, base::DataVector& alpha, double ylower,
                         bool limitNodalValues = true);
  void doUpperLimitation(base::Grid& grid, base::DataVector& alpha, double yupper,
                         bool limitNodalValues = true);
  void doLimitation(base::Grid& grid, base::DataVector& alpha, double ylower, double yupper,
                    bool limitNodalValues = true);

 private:
  void prepareForLowerLimitation(base::Grid& grid, base::DataVector& alpha, double ylower);
  void inverseFromLowerLimitation(base::Grid& grid, base::DataVector& alpha, double ylower);
  void prepareForUpperLimitation(base::Grid& grid, base::DataVector& alpha, double yupper);
  void inverseFromUpperLimitation(base::Grid& grid, base::DataVector& alpha, double yupper);

  void addConst(base::Grid& grid, base::DataVector& alpha, double c, double y);

  std::unique_ptr<datadriven::OperationMakePositive> opPositive;
  bool verbose;
};

} /* namespace datadriven */
} /* namespace sgpp */
