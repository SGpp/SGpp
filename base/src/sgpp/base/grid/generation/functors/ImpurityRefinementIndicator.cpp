// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/ImpurityRefinementIndicator.hpp>

//#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
//#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
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


ImpurityRefinementIndicator::ImpurityRefinementIndicator(Grid& grid, DataMatrix& dataset, 
  DataVector& alphas, DataVector& w1, DataVector& w2, DataVector& classesComputed,
  double threshold, size_t refinements_num) :
      grid(grid),
      dataset(dataset),
      alphas(alphas),    // for svm only
      w1(w1),            // for svm only
      w2(w2),            // for svm only
      classesComputed(classesComputed)
  {
    // set global variables accordingly
    this->refinementsNum = refinements_num;
    this->threshold = threshold;
}

double ImpurityRefinementIndicator::operator()(GridPoint& point) const {
  // the actual value of the impurity indicator
  double impurityInd = 0.0;
  // counter of datapoints lying on support of corresponding basis function
  size_t cnt = 0;

  size_t numClasses = 2; // ToDo: pass number of classes/labels and store in member
  DataVector fractions(numClasses);
  fractions.setAll(0);

  //SBasis& basis = const_cast<SBasis&>(grid.getBasis());

  // go through the whole dataset.
  for (size_t row = 0; row < dataset.getNrows(); ++row) {
    // level and index of the grid point in dimension d
    level_t level = 0;
    index_t index = 0;
    // mesh width in dimension d
    double h;

    double valueInDim;
    size_t cntSupDim = 0;

    for (size_t dim = 0; dim < point.getDimension(); ++dim) {
      level = point.getLevel(dim);
      index = point.getIndex(dim);
      h = std::pow(2, (-1.0)*static_cast<double>(level));

      valueInDim = dataset.get(row, dim);
      
      if ( (valueInDim >= (static_cast<double>(index) - 1.0)*h) && (valueInDim <= (static_cast<double>(index) + 1.0)*h) ) {
        cntSupDim ++;
      }
    }

    if (cntSupDim == point.getDimension()) {
      // determine current prediction (class) for datapoint
      double c = classesComputed.get(row);
      // increase respective fraction
      if (c == 1) {
        fractions.set(0, fractions.get(0)+1);
      }
      else {
        fractions.set(1, fractions.get(1)+1);
      }
      cnt ++;
    }
  }

  //impurityInd = 1.0;
  //impurityInd = 0.0;
  
  if (cnt > 0) {
    // Gini impurity
    for (size_t i = 0; i < numClasses; ++i) {
      impurityInd += std::pow(fractions.get(i)/static_cast<double>(cnt), 2); //relative
      //impurityInd += std::pow(fractions.get(i)/dataset.getNrows(), 2); //absolute
    }
    impurityInd = 1.0 - impurityInd;
  }

  return impurityInd;
}


/*double ImpurityRefinementIndicator::runOperator(GridStorage& storage,
    size_t seq) {
  return (*this)(storage.getPoint(seq));
}*/


size_t ImpurityRefinementIndicator::getRefinementsNum() const {
  return refinementsNum;
}

double ImpurityRefinementIndicator::getRefinementThreshold() const {
  return threshold;
}

double ImpurityRefinementIndicator::start() const {
  return 0.0;
}

double ImpurityRefinementIndicator::operator()(GridStorage& storage,
    size_t seq) const {
  throw std::logic_error("This form of the operator() is not implemented "
                         "for impurity indicators.");
}


void ImpurityRefinementIndicator::update(GridPoint& point) {
  SBasis& basis = const_cast<SBasis&>(grid.getBasis());

  double w1_new = 0.0;
  double w2_new = 0.0;
  double res;

  level_t level;
  index_t index;
  double value;
  double valueInDim;
  // go through whole dataset
  for (size_t row = 0; row < dataset.getNrows(); ++row) {
    // evaluate datapoint at current grid point in each dimension
    value = 1.0;
    for (size_t dim = 0; dim < point.getDimension(); ++dim) {
      level = point.getLevel(dim);
      index = point.getIndex(dim);
      valueInDim = dataset.get(row, dim);
      value *= basis.eval(level, index, valueInDim);
    }

    res = value * alphas.get(row); 
    w1_new += res;
    res = value * std::abs(alphas.get(row));
    w2_new += res;
  }

  w1.append(w1_new);
  w2.append(w2_new);
}



}  // namespace base
}  // namespace sgpp
