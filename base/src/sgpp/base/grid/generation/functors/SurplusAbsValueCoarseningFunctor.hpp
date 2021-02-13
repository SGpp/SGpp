// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * A refinement functor, refining according to the maximal absolute values in a DataVector provided,
 * weighted with the absolute value of the sparse grid approximation at the corresponding grid
 * point.
 */
class SurplusAbsValueCoarseningFunctor : public CoarseningFunctor {
 public:
  /**
   * Constructor.
   *
   * @param grid The grid on which we approximates the solution.
   * @param alpha DataVector that is basis for coarsening decisions. The i-th entry corresponds to
   * the i-th grid point.
   * @param removements_num Number of grid points which should be removed (if possible - there could
   * be less removable grid points), default: 1
   * @param threshold The absolute value of the entries have to be less or equal than the threshold
   * to be considered for coarsening, default: 0.0
   */
  SurplusAbsValueCoarseningFunctor(Grid& grid, DataVector& alpha, size_t removements_num = 1,
                                   double threshold = 0.0);

  /**
   * Destructor
   */
  ~SurplusAbsValueCoarseningFunctor() override;

  double operator()(GridStorage& storage, size_t seq) override;

  double start() const override;

  size_t getRemovementsNum() const override;

  double getCoarseningThreshold() const override;

 protected:
  /// pointer to the grid
  Grid& grid;

  /// pointer to the vector that stores the alpha values
  DataVector& alpha;

  /// evaluation of the grid at gridpoints
  DataVector evals;

  /// number of grid points to remove
  size_t removements_num;

  /**
   * threshold, only the points with greater to equal absolute values of the
   * refinement criterion (e.g. alpha or error) will be refined
   */
  double threshold;
};

}  // namespace base
}  // namespace sgpp
