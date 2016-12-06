// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef IMPURITYREFINEMENTINDICATOR_HPP_
#define IMPURITYREFINEMENTINDICATOR_HPP_

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>

#include <sgpp/globaldef.hpp>

#include <utility>

namespace sgpp {
namespace base {

/**
 *  A refinement indicator for classification problems based on impurity
 * measures
 *  (e.g. gini impurity, entropy impurity,...). It calculates local impurities
 * based
 *  on the information from the provided data set. If the indicator is applied
 * within
 *  the SVM learner, the normal vector needs to be extended after each
 * refinement.
 */

class ImpurityRefinementIndicator : public RefinementFunctor {
 public:
  typedef std::pair<size_t, double> value_type;

  typedef GridPoint counter_key_type;

  /**
   * Constructor.
   *
   * @param grid The grid to refine.
   * @param dataset The set of data points used to compute impurities
   * @param alphas The weights corresponding to the support vectors (only
   * required for SVM learner)
   * @param w1 Normal vector (only required for SVM learner)
   * @param w2 Normal vector computed with abs values (only required for SVM
   * learner)
   * @param classesComputed The predicted labels for the data points from
   * dataset
   * @param threshold The refinement threshold; Only grid points with
   *        indicator values greater than this threshold will be refined
   * @param refinementsNum The max amount of grid points to be refined
   */
  ImpurityRefinementIndicator(Grid& grid, DataMatrix& dataset,
                              DataVector* alphas, DataVector* w1,
                              DataVector* w2, DataVector& classesComputed,
                              double threshold = 0.0,
                              size_t refinementsNum = 1);

  /**
   * Destructor
   */
  virtual ~ImpurityRefinementIndicator() {}

  /**
   * This should be returning a refinement indicator for the specified grid
   * point.
   * The point with the highest value will be refined first.
   *
   * @param point The grid point for which to calculate an indicator value
   * @return The indicator value
   */
  virtual double operator()(GridPoint& point) const;

  /**
   * Returns the maximal number of points that should be refined.
   * The maximal number of points to refine is set in the constructor of the
   * implementing class.
   *
   * @return Number of points that should refined. Default value: 1.
   */
  size_t getRefinementsNum() const override;

  /**
   * Returns the threshold for refinement.
   * Only the grid points with absolute value of refinement criterion greater
   * than this threshold will be refined.
   *
   * @return Threshold value for refinement. Default value: 0.
   */
  double getRefinementThreshold() const override;

  double start() const override;

  /**
  * This should be returning a refinement value for every grid point.
  * The point with the highest value will be refined first.
  *
  * @param storage Reference to the grids storage object
  * @param seq Sequence number in the coefficients array
  *
  * @return refinement value
  */
  double operator()(GridStorage& storage, size_t seq) const override;

  /**
   * Update normal vector of SVM. For each new grid point the normal vector
   * has to be extended by one component. Only required for SVMLearner!
   *
   * @param point The new grid point
   */
  void update(GridPoint& point);

  DataVector* alphas;  // required for svm learner only

  DataVector* w1;  // required for svm learner only

  DataVector* w2;  // required for svm learner only

 protected:
  // data set that will be evaluated
  DataMatrix& dataset;
  // the corresponding computed class labels
  DataVector& classesComputed;
  // max number of grid points to refine
  size_t refinementsNum;
  // threshold, only the points with greater absolute values
  // of the refinement criterion will be refined
  double threshold;

 private:
  // the grid to refine
  Grid& grid;
};

}  // namespace base
}  // namespace sgpp

#endif /* IMPURITYREFINEMENTINDICATOR_HPP_ */
