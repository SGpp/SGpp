// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/SqXdPhidPhiUpBBLinearStretched.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

SqXdPhidPhiUpBBLinearStretched::SqXdPhidPhiUpBBLinearStretched(sgpp::base::GridStorage* storage)
    : storage(storage), stretching(storage->getStretching()) {}

SqXdPhidPhiUpBBLinearStretched::~SqXdPhidPhiUpBBLinearStretched() {}

void SqXdPhidPhiUpBBLinearStretched::operator()(sgpp::base::DataVector& source,
                                                sgpp::base::DataVector& result,
                                                grid_iterator& index, size_t dim) {
  // get boundary values
  double fl = 0.0;
  double fr = 0.0;

  rec(source, result, index, dim, fl, fr);
}

void SqXdPhidPhiUpBBLinearStretched::rec(sgpp::base::DataVector& source,
                                         sgpp::base::DataVector& result, grid_iterator& index,
                                         size_t dim, double& fl, double& fr) {
  size_t seq = index.seq();

  fl = fr = 0.0;
  double fml = 0.0;
  double fmr = 0.0;

  sgpp::base::GridStorage::index_type::level_type current_level;
  sgpp::base::GridStorage::index_type::index_type current_index;

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
  double posl = 0, posr = 0, posc = 0;
  this->stretching->getAdjacentPositions(static_cast<int>(current_level),
                                         static_cast<int>(current_index), dim, posc, posl, posr);
  double baseLength = posr - posl;
  double leftLength = posc - posl;
  double rightLength = posr - posc;

  double fm = fml + fmr;

  double alpha_value = source[seq];

  double c = 1.0 / 3.0 * (posc + posr + posl);
  // transposed operations:
  result[seq] = fm;

  fl = c * alpha_value + fl + fm * (rightLength / baseLength);
  fr = -c * alpha_value + fr + fm * (leftLength / baseLength);
}

}  // namespace finance
}  // namespace sgpp
