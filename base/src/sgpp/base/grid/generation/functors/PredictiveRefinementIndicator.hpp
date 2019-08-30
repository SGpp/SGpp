// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PREDICTIVEREFINEMENTINDICATOR_HPP_
#define PREDICTIVEREFINEMENTINDICATOR_HPP_


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
 *  A refinement error indicator for regression problems based on the residuals of the datasets.
 *
 *  It calculates an error messure based on the information from the data set:
 *  For a new grid point g on level l and index i, it calculates the indicator as a
 *  sum of squared residuals ( (value of sample from dataset - grid evaluated at the sample's coordinates), squared),
 *  weighted with the underlying basis function of the grid point.
 */

class PredictiveRefinementIndicator: public RefinementFunctor {
 public:
  typedef std::pair<size_t, double> value_type;

  typedef GridPoint counter_key_type;
  typedef uint64_t counter_value_type;

  /**
   * Constructor.
   *
   * @param grid DataVector that is basis for refinement decisions. The i-th entry corresponds to the i-th grid point.
   * @param dataSet contains all points of the source data set. Each row contains coordinates of a single grid point,
   * without the function evaluation (meaning only data from omega).
   * @param errorVector a DataVector containing the squared absolute error
   * (given value of data point - evaluation of sparse grid at the data point position) for each grid point in dataSet.
   * @param refinements_num the amount of grid points to maximally be refined or created, depending on refinement strategy.
   * @param threshold The absolute value of the entries have to be greater or equal than the threshold
   * @param minSupportPoints The minimal number of data points that have to be within the support of a basis function
   * for refinement.
   */
  PredictiveRefinementIndicator(Grid& grid, DataMatrix& dataSet,
                                DataVector& errorVector,
                                size_t refinements_num = 1,
                                double threshold = 0.0,
                                uint64_t minSupportPoints = 0);


  /**
   * Destructor
   */
  ~PredictiveRefinementIndicator() override {}


  /**
   * This should be returning a refinement indicator for the specified grid point
   * The point with the highest value will be refined first.
   *
   * @param point for which to calculate an indicator value
   * @return refinement value
   */
  virtual double operator()(GridPoint& point) const;

  double runOperator(GridStorage& storage, size_t seq);


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
   * Only the grid points with absolute value of refinement criterion (e.g., alpha or error) greater
   * or equal to this threshold will be refined.
   *
   * @return threshold value for refinement. Default value: 0.
   */
  double getRefinementThreshold() const override;


  double start() const override;

  /**
   * Returns the lower bound of refinement criterion (e.g., alpha or error
   */
  uint64_t getMinSupportPoints() const {
    return minSupportPoints_;
  }

  void setMinSupportPoints(uint64_t minSupportPoints) {
    minSupportPoints_ = minSupportPoints;
  }

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

 protected:
  // for each Point in the dataSet, this Contains the squared absolute
  // offset between sparse grid and data point value;
  DataVector& errorVector;
  /// number of grid points to refine
  size_t refinementsNum;
  /// threshold, only the points with greater to equal absolute values
  // of the refinement criterion (e.g. alpha or error) will be refined
  double threshold;

  // data set that will be evaluated
  DataMatrix& dataSet;

  /*
   * Evaluates a basis function at a point on level "level" and index "index". The type of basis function is determinded
   * by the grid type, who's integer representation is saved in the member gridType.
   * @param level on which the basis function is located
   * @param index within the level "level"
   * @param value represents the x value on which the grid should be evaluated
   * @return basis function of grid type evaluated at "value" on level "level" and index "index".
   */
  double basisFunctionEvalHelper(level_t level, index_t index, double value);

 private:
  /*
   * Each grid point has a support, that is constructed via a tensor product approach.
   * This function is a helper method for the operator()(index_type*), determining min(supp(basisFunction(level,index))) and max(supp(BasisFunction(level,index))) of the basis function associated with the grid type
   * @param point for which to calculate the support Vector
   * @param dataVector in the dimensions of the grid Point, where each row of the data vector represents the minimum of the
   * support in the dimension of the grid point
   * @param dataVector in the dimensions of the grid Point, where each row of the data vector represents the maximum of the
   * support in the dimension of the grid point.
   */
  void buildGPSupportMask(GridPoint& point, DataVector& floorMask, DataVector& ceilingMask);

  /*
   * Figures out if an data point from the data set is on the support of the basis function associated with a grid point. Helper method for the operator()(index_type*)
   * @param floorMask DataVector which holds min(supp(basisFunction)) for one basis functions
   * @param ceilingMask DataVector which holds max(supp(BasisFunction)) for one basis function
   * @param entry size_t row in data set to be analyzed
   * @return true if the point from the dataset is located on the support of the basis function, else false.
   */
  bool isOnSupport(DataVector& floorMask, DataVector& ceilingMask,
                   size_t entry);

  /*
   * integer representation of the grid type needed for evaluation of basis functions.
   */
  sgpp::base::GridType gridType;

  Grid& grid_;

  uint64_t minSupportPoints_;
};

}  // namespace base
}  // namespace sgpp
#endif /* PREDICTIVEREFINEMENTINDICATOR_HPP_ */
