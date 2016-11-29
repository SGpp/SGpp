// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef FORWARDSELECTORREFINEMENTINDICATOR_HPP_
#define FORWARDSELECTORREFINEMENTINDICATOR_HPP_

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace base {

/**
 *  A refinement indicator for support vector classification using 
 *  sparse grids (according to König BA). 
 */

class ForwardSelectorRefinementIndicator: public RefinementFunctor {
 public:
  typedef std::pair<size_t, double> value_type;

  typedef GridPoint counter_key_type;

  /**
   * Constructor.
   *
   * @param grid The sparse grid
   * @param svs Contains all currently stored support vectors
   * @param alphas The weights corresponding to the support vectors
   * @param w1 The normal vector
   * @param w2 The normal vector computted with abs weights
   * @param beta Specifies relevance of grid points (default: equal relevance for all grid points)
   * @param threshold The refinement threshold; Only grid points with 
   *        indicator values greater than this threshold will be refined
   * @param refinementsNum The max amount of grid points to be refined  
   */
  ForwardSelectorRefinementIndicator(Grid& grid, DataMatrix& svs,
                                     DataVector& alphas,
                                     DataVector& w1,
                                     DataVector& w2,
                                     double beta,
                                     double threshold = 0.0,
                                     size_t refinementsNum = 1);

  /**
   * Destructor
   */
  virtual ~ForwardSelectorRefinementIndicator() {}
  
  /**
  * This should be returning a refinement indication value for every grid point.
  * The grid point with the highest value will be refined first.
  *
  * @param storage Reference to the grids storage object
  * @param seq Sequence number in the coefficients array
  * @return The refinement indicator value
  */
  double operator()(GridStorage& storage, size_t seq) const override;

  double runOperator(GridStorage& storage, size_t seq);

  /**
   * Returns the maximal number of points that should be refined.
   * The maximal number of points to refine is set in the constructor of the implementing class.
   *
   * @return number of points that should refined. Default value: 1.
   */
  size_t getRefinementsNum() const override;

  /**
   * Returns the threshold for refinement.
   * Only the grid points with absolute value of refinement criterion greater
   * than this threshold will be refined.
   *
   * @return threshold value for refinement. Default value: 0.
   */
  double getRefinementThreshold() const override;

  double start() const override;

  /**
   * This should be returning a refinement indicator for the specified grid point
   * The point with the highest value will be refined first.
   *
   * @param point Grid point for which to calculate an indicator value
   * @return The indicator value
   */
  virtual double operator()(GridPoint& point) const;

  /**
   * Update normal vector of SVM. For each new grid point the normal vector
   * has to be extended by one component.
   *
   * @param point The new grid point 
   */
  void update(GridPoint& point);

 protected:
  // set of support vectors that will be evaluated
  DataMatrix& svs;
  // the normal vector
  DataVector& w1;
  // the normal vector computed with absolute weights
  DataVector& w2;
  // the weights corresponding to the support vectors
  DataVector& alphas;
  // current loss
  std::shared_ptr<DataVector> rv1; 
  // current loss (abs)
  std::shared_ptr<DataVector> rv2;
  // specifies relevance of grid points
  double beta;
  // max number of grid points to refine
  size_t refinementsNum;
  // threshold, only the points with greater indicator values be refined
  double threshold;

 private:
  // the sparse grid
  Grid& grid;

};

}  // namespace base
}  // namespace sgpp

#endif /* FORWARDSELECTORREFINEMENTINDICATOR_HPP_ */
