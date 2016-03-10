// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

PhiPhiUpBBLinearStretched::PhiPhiUpBBLinearStretched(sgpp::base::GridStorage* storage)
    : storage(storage), stretching(storage->getStretching()) {}

PhiPhiUpBBLinearStretched::~PhiPhiUpBBLinearStretched() {}

void PhiPhiUpBBLinearStretched::operator()(sgpp::base::DataVector& source,
                                           sgpp::base::DataVector& result, grid_iterator& index,
                                           size_t dim) {
  // get boundary values
  double fl = 0.0;
  double fr = 0.0;

  rec(source, result, index, dim, fl, fr);
}

void PhiPhiUpBBLinearStretched::rec(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                    grid_iterator& index, size_t dim, double& fl, double& fr) {
  size_t seq = index.seq();

  fl = fr = 0.0;
  double fml = 0.0;
  double fmr = 0.0;

  sgpp::base::GridStorage::index_type::level_type current_level;
  sgpp::base::GridStorage::index_type::index_type current_index;

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->end(index.seq())) {
      rec(source, result, index, dim, fl, fml);
    }

    index.stepRight(dim);

    if (!storage->end(index.seq())) {
      rec(source, result, index, dim, fmr, fr);
    }

    index.up(dim);
  }

  index.get(dim, current_level, current_index);
  double posl = 0, posr = 0, currentPosition = 0;
  this->stretching->getAdjacentPositions(static_cast<int>(current_level),
                                         static_cast<int>(current_index), dim, currentPosition,
                                         posl, posr);
  double baseLength = posr - posl;
  double leftLength = currentPosition - posl;
  double rightLength = posr - currentPosition;

  double fm = fml + fmr;

  double alpha_value = source[seq];

  // transposed operations:
  result[seq] = fm;

  fl =
      (1.0 / 6.0) * (baseLength + rightLength) * alpha_value + fl + fm * (rightLength / baseLength);
  fr = (1.0 / 6.0) * (baseLength + leftLength) * alpha_value + fr + fm * (leftLength / baseLength);
}

}  // namespace pde
}  // namespace sgpp
