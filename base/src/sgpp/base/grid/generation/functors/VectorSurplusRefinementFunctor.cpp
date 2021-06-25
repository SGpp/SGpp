// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/VectorSurplusRefinementFunctor.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

VectorSurplusRefinementFunctor::VectorSurplusRefinementFunctor(DataMatrix& alphas,
                                                               size_t refinements_num,
                                                               double threshold)
    : alphas(alphas), refinements_num(refinements_num), threshold(threshold) {}

VectorSurplusRefinementFunctor::~VectorSurplusRefinementFunctor() {}

double VectorSurplusRefinementFunctor::operator()(GridStorage& storage, size_t seq) const {
  sgpp::base::DataVector row(alphas.getNcols());
  alphas.getRow(seq, row);
  row.abs();
  double val = row.max();
  return val;
}

double VectorSurplusRefinementFunctor::start() const { return 0.0; }

size_t VectorSurplusRefinementFunctor::getRefinementsNum() const { return this->refinements_num; }

double VectorSurplusRefinementFunctor::getRefinementThreshold() const { return this->threshold; }

}  // namespace base
}  // namespace sgpp
