// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <map>
#include <string>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <algorithm>


namespace sgpp {
namespace base {


PredictiveRefinementIndicator::PredictiveRefinementIndicator(Grid& grid,
    DataMatrix& dataSet, DataVector& errorVector,
    size_t refinements_num, double threshold,
    uint64_t minSupportPoints):
        errorVector(errorVector), dataSet(dataSet), grid_(grid),
        minSupportPoints_(minSupportPoints) {
  // find out what type of grid is used;
  gridType = grid.getType();

  // set global Variables accordingly
  this->refinementsNum = refinements_num;
  this->threshold = threshold;
}

double PredictiveRefinementIndicator::operator ()(
  AbstractRefinement::index_type& gridPoint) const {
  // the actual value of the errorIndicator
  double errorIndicator = 0.0;
  double denominator = 0.0;
  double r22 = 0.0;
  double r2phi = 0.0;


  // counter of contributions - for DEBUG purposes
  size_t counter = 0;

  SBasis& basis = const_cast<SBasis&>(grid_.getBasis());
  // go through the whole dataset.
  // -> if data point on the support of the grid point in all
  // dim then calculate error Indicator.
  #pragma omp parallel for schedule(static) \
  reduction(+:errorIndicator, denominator, r22, r2phi, counter)

  for (size_t row = 0; row < dataSet.getNrows(); ++row) {
    // level, index and evaulation of a gridPoint in dimension d
    AbstractRefinement::level_t level = 0;
    AbstractRefinement::index_t index = 0;
    double valueInDim;
    // if on the support of the grid point in all dim
    // if(isOnSupport(&floorMask,&ceilingMask,row))
    // {
    // counter for DEBUG
    // ++counter;*****
    double funcval = 1.0;

    // calculate error Indicator
    for (size_t dim = 0; dim < gridPoint.getDimension() && funcval != 0; ++dim) {
      level = gridPoint.getLevel(dim);
      index = gridPoint.getIndex(dim);

      valueInDim = dataSet.get(row, dim);

      funcval *=  std::max(0.0, basis.eval(level, index, valueInDim));

      // basisFunctionEvalHelper(level,index,valueInDim);
    }

    errorIndicator += funcval * errorVector.get(row)/**errorVector->get(row)*/;
    r22 += errorVector.get(row) * errorVector.get(row);
    r2phi += funcval * errorVector.get(row);
    denominator += funcval * funcval;

    if (funcval != 0.0) counter++;

    // }
  }

  AbstractRefinement::index_type idx(gridPoint);

  if (denominator != 0 && counter >= minSupportPoints_) {
    // to match with OnlineRefDim, use this:
    // return (errorIndicator * errorIndicator) / denominator;

    double a = (errorIndicator / denominator);
    return /*r2phi/denominator*/ /*2*r22 - 2*a*r2phi + a*a*denominator*/ a *
        (2 * r2phi - a * denominator);
    // return fabs(a);
  } else {
    return 0.0;
  }
}


/*double PredictiveRefinementIndicator::operator ()(GridStorage& storage, size_t seq) {
  return errorVector->get(seq);
}*/


double PredictiveRefinementIndicator::runOperator(GridStorage& storage,
    size_t seq) {
  return (*this)(storage.getGridIndex(seq));
}


size_t PredictiveRefinementIndicator::getRefinementsNum() const {
  return refinementsNum;
}

double PredictiveRefinementIndicator::getRefinementThreshold() const {
  return threshold;
}

double PredictiveRefinementIndicator::start() const {
  return 0.0;
}

double PredictiveRefinementIndicator::operator()(GridStorage& storage,
    size_t seq) const {
  throw std::logic_error("This form of the operator() is not implemented "
                         "for predictive indicators.");
}




}  // namespace base
}  // namespace sgpp
