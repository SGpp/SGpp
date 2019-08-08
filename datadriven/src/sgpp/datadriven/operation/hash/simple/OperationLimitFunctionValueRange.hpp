// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationMakePositive.hpp>

#include <vector>
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
      bool verbose = false, sgpp::base::ScalarFunction* f = nullptr);

  virtual ~OperationLimitFunctionValueRange();

  void doLowerLimitation(base::Grid& grid, base::DataVector& alpha, double ylower,
                         bool limitNodalValues = true);
  void doUpperLimitation(base::Grid& grid, base::DataVector& alpha, double yupper,
                         bool limitNodalValues = true);
  void doLimitation(base::Grid& grid, base::DataVector& alpha, double ylower, double yupper,
                    bool limitNodalValues = true);

  /**
   *
   * @return vector containing the indices of the added grid points
   */
  std::vector<size_t>& getAddedGridPoints();

  /**
   *
   * @return vector containing the indices which have just been added for range limiting
   */
  std::vector<size_t>& getAddedGridPointsForRangeLimitation();

  /**
   *
   * @return number of newly added grid points
   */
  size_t numAddedGridPoints();

  /**
   *
   * @return number of newly added grid points for guaranteeing the range
   */
  size_t numAddedGridPointsForRangeLimitation();

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
