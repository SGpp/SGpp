// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationMakePositiveCandidateSetAlgorithm.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationMakePositiveInterpolationAlgorithm.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>
#include <map>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>

namespace sgpp {
namespace datadriven {

enum class MakePositiveCandidateSearchAlgorithm {
  FullGrid,
  Intersections,
  HybridFullIntersections,
  IntersectionsJoin
};
enum class MakePositiveInterpolationAlgorithm {
  SetToZero,
  InterpolateExp,
  InterpolateBoundaries1d,
  InterpolateFunction
};

/**
 * This class enforces the function value range of a sparse grid function to be larger than 0.
 * It uses a discretization based approach where we add the minimum amount of full grid points
 * to enforce the positivity.
 */
class OperationMakePositive {
 public:
  typedef std::map<size_t, base::HashGridPoint> gridPointCandidatesMap;

  /**
   * Constructor.
   *
   * @param candidateSearchAlgorithm defines how to generate the full grid candidate set
   * @param interpolationAlgorithm defines how to compute the coefficients for the new grid points
   * @param generateConsistentGrid define if the hierarchical ancestors of all new grid points are
   * inserted as well
   * @param verbose print information or not
   * @param f scalar function to be interpolated
   */
  explicit OperationMakePositive(MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm =
                                     MakePositiveCandidateSearchAlgorithm::IntersectionsJoin,
                                 MakePositiveInterpolationAlgorithm interpolationAlgorithm =
                                     MakePositiveInterpolationAlgorithm::SetToZero,

                                 bool generateConsistentGrid = true, bool verbose = false,
                                 sgpp::base::ScalarFunction* f = nullptr);

  /**
   * Descrutor
   */
  virtual ~OperationMakePositive();

  /**
   * initializes the operation
   *
   * @param grid Grid
   * @param alpha coefficients
   */
  void initialize(base::Grid& grid, base::DataVector& alpha);

  /**
   * Make the sparse grid function defined by grid and coefficient vector positive.
   *
   * @param grid Grid where the new grid is stored
   * @param alpha coefficient vector of new grid
   * @param forcePositiveNodalValues nodal values are forced to be non-negative or not
   */
  void makePositive(base::Grid& grid, base::DataVector& alpha,
                    bool forcePositiveNodalValues = false);

  /**
   * Enforce the function values at each grid point to larger than the specified tolerance. The ones
   * which are not are set to zero. For this function we need the hierarchization and
   * dechierarchization operations.
   *
   * @param grid grid
   * @param alpha coefficient vector
   * @param tol tolerance for negativity
   */
  void makeCurrentNodalValuesPositive(base::Grid& grid, base::DataVector& alpha,
                                      double tol = -1e-14);

  /**
   *
   * @return vector containing the indices of the added grid points
   */
  std::vector<size_t>& getAddedGridPoints();

  /**
   *
   * @return vector containing the indices which have just been added for positivity
   */
  std::vector<size_t>& getAddedGridPointsForPositivity();

  /**
   *
   * @return number of newly added grid points
   */
  size_t numAddedGridPoints();

  /**
   *
   * @return number of newly added grid points for guaranteeing positivity
   */
  size_t numAddedGridPointsForPositivity();

  /**
   *
   * @return number of newly added grid points for guaranteeing
   * positivity per iteration
   */
  base::DataVector& numAddedGridPointsForPositivityPerIteration();

  /**
   * @return candidate set algorithm
   */
  OperationMakePositiveCandidateSetAlgorithm& getCandidateSetAlgorithm();

 private:
  /**
   * Enforce the function values at the new grid points to be positive. It is similar to the
   * ´makeCurrentNodalValuesPositive´ but does the hierarchization directly.
   *
   * @param grid Grid
   * @param alpha coefficient vector
   * @param newGridPoints new grid points for which we need to compute the coefficients
   * @param tol tolerance for positivity
   */
  void forceNewNodalValuesToBePositive(base::Grid& grid, base::DataVector& alpha,
                                       std::vector<size_t>& newGridPoints, double tol = -1e-14);

  /**
   * Extract the non existing candidates from the candidate set which have the currently checked
   * level sum.
   *
   * @param newGrid Grid
   * @param candidates candidate set
   * @param currentLevelSum currently checked grid points with this level sum
   * @param finalCandidates candidate set where all the grid points do not exist in newGrid and
   * their level sum is equal to currentLevelSum
   */
  void extractNonExistingCandidatesByLevelSum(
      base::Grid& newGrid, std::vector<std::shared_ptr<base::HashGridPoint>>& candidates,
      size_t currentLevelSum, std::vector<std::shared_ptr<base::HashGridPoint>>& finalCandidates);

  /**
   * Adds all the candidates from which the function value is smaller than the tolerance.
   *
   * @param grid Grid
   * @param alpha coefficient vector
   * @param candidates candidate set
   * @param addedGridPoints newly added grid points
   * @param tol tolerance
   */
  void addFullGridPoints(base::Grid& grid, base::DataVector& alpha,
                         std::vector<std::shared_ptr<base::HashGridPoint>>& candidates,
                         std::vector<size_t>& addedGridPoints, double tol = -1e-14);

  /// maximum level of the reference full grid
  size_t maxLevel;

  /// range for level sums to be tested
  size_t minimumLevelSum;
  size_t maximumLevelSum;

  /// needed new grid points for positivity
  std::vector<size_t> addedGridPointsForPositivity;
  /// newly added grid points
  std::vector<size_t> addedGridPoints;
  /// counter for the number of added grid points for positivity in each iteration
  base::DataVector countAddedGridPointsForPositivityPerIteration;
  /// sets if a consistent grid is computed or not
  bool generateConsistentGrid;

  /// candidate search algorithm
  datadriven::MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm;
  datadriven::MakePositiveInterpolationAlgorithm interpolationAlgorithm;

  std::shared_ptr<datadriven::OperationMakePositiveCandidateSetAlgorithm> candidateSearch;
  std::shared_ptr<datadriven::OperationMakePositiveInterpolationAlgorithm> interpolation;
  sgpp::base::ScalarFunction* f;

  /// verbosity
  bool verbose;
};

} /* namespace datadriven */
} /* namespace sgpp */
