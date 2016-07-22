// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef IMPURITYREFINEMENTINDICATOR_HPP_
#define IMPURITYREFINEMENTINDICATOR_HPP_


#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>


#include <sgpp/globaldef.hpp>


#include <unordered_map>
#include <utility>


namespace sgpp {
namespace base {

/**
 *  A refinement indicator for classification problems based on impurity measures
 *  (gini impurity, entropy impurity etc.).
 *  It calculates local impurities based on the information from the dataset:
 *  For a grid point g on level l and index i, it calculates the indicator as a
 *  
 * 
 */

class ImpurityRefinementIndicator: public RefinementFunctor {
 public:
  typedef std::pair<size_t, double> value_type;

  typedef GridPoint counter_key_type;

  /**
   * Constructor.
   *
   * @param grid DataVector that is basis for refinement decisions. The i-th entry corresponds to the i-th grid point.
   * @param dataset contains all points of the source data set.
   */
  ImpurityRefinementIndicator(Grid& grid, DataMatrix& dataset,
                              DataVector& alphas, DataVector& w1, DataVector& w2,
                              DataVector& classesComputed,
                              double threshold = 0.0,
                              size_t refinements_num = 1);


  /**
   * Destructor
   */
  virtual ~ImpurityRefinementIndicator() {}


  /**
   * This should be returning a refinement indicator for the specified grid point
   * The point with the highest value will be refined first.
   *
   * @param point for which to calculate an indicator value
   * @return refinement value
   */
  virtual double operator()(GridPoint& point) const;

  //double runOperator(GridStorage& storage, size_t seq);


  /**
   * Returns the maximal number of points that should be refined.
   *
   * The maximal number of points to refine is set in the constructor of the implementing class.
   *
   * @return number of points that should refined. Default value: 1.
   */
  size_t getRefinementsNum() const override;

  /**
   * Returns the threshold for refinement.
   *
   * Only the grid points with absolute value of refinement criterion greater
   * or equal to this threshold will be refined.
   *
   * @return threshold value for refinement. Default value: 0.
   */
  double getRefinementThreshold() const override;


  double start() const override;

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

  void update(GridPoint& point);

 protected:
  // data set that will be evaluated
  DataMatrix& dataset;

  DataVector& alphas; // for svm only

  DataVector& w1;     // for svm only

  DataVector& w2;     // for svm only

  DataVector& classesComputed;

  /// number of grid points to refine
  size_t refinementsNum;

  /// threshold, only the points with greater to equal absolute values
  // of the refinement criterion will be refined
  double threshold;

 private:
  /*
   * integer representation of the grid type needed for evaluation of basis functions.
   */
  //sgpp::base::GridType gridType;

  Grid& grid;
};

}  // namespace base
}  // namespace sgpp
#endif /* IMPURITYREFINEMENTINDICATOR_HPP_ */
