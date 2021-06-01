// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * A refinement functor for vector valued function.
 * Refining according to the maximal absolute values in a DataMatrix (of coefficients) provided.
 */
class VectorSurplusRefinementFunctor : public RefinementFunctor {
 public:
  /**
   * Constructor.
   *
   * @param alphas DataMatrix that is basis for refinement decisions. The [i,j]-th entry
   corresponds
   * to the i-th grid point and the j-th result (<=> each column are coefficients for one output
   * dimension)
   * @param refinements_num Number of grid points which should be refined (if possible - there
   could
   * be less refinable grid points)
   * @param threshold The absolute value of the entries have to be greater or equal than the
   * threshold
   */
  VectorSurplusRefinementFunctor(DataMatrix& alphas, size_t refinements_num = 1,
                                 double threshold = 0.0);

  /**
   * Destructor
   */
  ~VectorSurplusRefinementFunctor() override;

  /**
   * This should be returning a refinement value for every grid point.
   * The point with the highest value will be refined first.
   *
   * @param storage reference to the grids storage object
   * @param seq sequence number in the coefficients array
   *
   * @return refinement value
   */
  double operator()(GridStorage& storage, size_t seq) const override;

  /**
   * Returns the lower bound of refinement criterion (e.g., alpha or error) (lower bound).
   * The refinement value of grid points to be refined have to be larger than this value
   *
   * @return lower bound
   */
  double start() const override;

  size_t getRefinementsNum() const override;

  double getRefinementThreshold() const override;

 protected:
  /// pointer to the vector that stores the alpha values
  DataMatrix& alphas;

  /// number of grid points to refine
  size_t refinements_num;

  /**
   * threshold, only the points with greater to equal absolute values of the
   * refinement criterion (e.g. alpha or error) will be refined
   */
  double threshold;
};

}  // namespace base
}  // namespace sgpp
