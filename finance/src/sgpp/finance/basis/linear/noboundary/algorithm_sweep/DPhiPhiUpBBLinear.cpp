// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

DPhiPhiUpBBLinear::DPhiPhiUpBBLinear(sgpp::base::GridStorage* storage)
    : storage(storage), boundingBox(storage->getBoundingBox()) {}

DPhiPhiUpBBLinear::~DPhiPhiUpBBLinear() {}

void DPhiPhiUpBBLinear::operator()(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                   grid_iterator& index, size_t dim) {
  // get boundary values
  double fl = 0.0;
  double fr = 0.0;

  rec(source, result, index, dim, fl, fr);
}

void DPhiPhiUpBBLinear::rec(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                            grid_iterator& index, size_t dim, double& fl, double& fr) {
  size_t seq = index.seq();

  fl = fr = 0.0;
  double fml = 0.0;
  double fmr = 0.0;

  sgpp::base::level_t current_level;
  sgpp::base::index_t current_index;

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fl, fml);
    }

    index.stepRight(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fmr, fr);
    }

    index.up(dim);
  }

  index.get(dim, current_level, current_index);

  double fm = fml + fmr;

  double alpha_value = source[seq];

  // transposed operations:
  result[seq] = fm;

  fl = (fm / 2.0) + (alpha_value * 0.5 + fl);
  fr = (fm / 2.0) + (alpha_value * (-0.5) + fr);
}

}  // namespace finance
}  // namespace sgpp
