// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/ForwardSelectorRefinementIndicator.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>

#include <map>
#include <string>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <algorithm>


namespace sgpp {
namespace base {


ForwardSelectorRefinementIndicator::ForwardSelectorRefinementIndicator(Grid& grid, DataMatrix& svs, 
    DataVector& alphas, DataVector& w1, DataVector& w2, double beta, double threshold, size_t refinements_num) :
        grid(grid),
        svs(svs),
        alphas(alphas),
        w1(w1),
        w2(w2),
        beta(beta),
        rv1(new DataVector(grid.getSize(), 0.0)),
        rv2(new DataVector(grid.getSize(), 0.0)) 
 {
  // set global Variables accordingly
  this->refinementsNum = refinements_num;
  this->threshold = threshold;
  
  // compute current loss using set of support vectors  
  DataVector losses(svs.getNrows()); // define as member of indicator if needed again later
  for (size_t i = 0; i < svs.getNrows(); i++) {
    DataVector x(svs.getNcols());
    svs.getRow((size_t)i, x);
    // compute transformation of current sv -> maybe call opSGKernel here ?
    DataVector xTrans(grid.getSize()); 
    DataMatrix xMatrix(1,svs.getNcols());
    xMatrix.setRow(0,x);
    DataVector alpha(1,1.0);
    op_factory::createOperationMultipleEval(grid, xMatrix)->multTranspose(alpha, xTrans); 
    double t = alphas.get(i);
    losses.set(i, 1.0 - w1.dotProduct(xTrans) * t);

    if (losses.get(i) > 0) {
      rv2->add(xTrans);
      xTrans.mult(t);
      rv1->add(xTrans);
    }
  }

}

double ForwardSelectorRefinementIndicator::operator()(GridStorage& storage, size_t seq) const {
  double epsilon = 0.000001;
  double measure = std::abs(w1.get(seq) + epsilon) / (std::pow(w2.get(seq),beta) * rv2->get(seq) + epsilon); 
  
  //std::cout << "measure: " << measure << std::endl;

  return measure;
}

double ForwardSelectorRefinementIndicator::runOperator(GridStorage& storage,
    size_t seq) {
  return (*this)(storage.getPoint(seq));
}

size_t ForwardSelectorRefinementIndicator::getRefinementsNum() const {
  return refinementsNum;
}

double ForwardSelectorRefinementIndicator::getRefinementThreshold() const {
  return threshold;
}

double ForwardSelectorRefinementIndicator::start() const {
  return 0.0;
}

double ForwardSelectorRefinementIndicator::operator()(GridPoint& point) const {
  throw std::logic_error("This form of the operator() is not implemented "
                         "for svm indicators.");
}

//void ForwardSelectorRefinementIndicator::update(ForwardSelectorRefinement::insertables_container_type& insertables) {
void ForwardSelectorRefinementIndicator::update(GridPoint& point) {
  SBasis& basis = const_cast<SBasis&>(grid.getBasis());

  double w1_new = 0.0;
  double w2_new = 0.0;
  double res;

  level_t level;
  index_t index;
  double value;
  double valueInDim;
  // go through all support vectors
  for (size_t row = 0; row < svs.getNrows(); ++row) {
    // evaluate sv at current grid point in each dimension
    value = 1.0;
    for (size_t dim = 0; dim < point.getDimension(); ++dim) {
      level = point.getLevel(dim);
      index = point.getIndex(dim);
      valueInDim = svs.get(row, dim);
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
