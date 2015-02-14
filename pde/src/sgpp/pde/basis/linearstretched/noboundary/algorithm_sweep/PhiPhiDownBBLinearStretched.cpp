// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {



    PhiPhiDownBBLinearStretched::PhiPhiDownBBLinearStretched(SGPP::base::GridStorage* storage) : storage(storage), stretching(storage->getStretching()) {
    }

    PhiPhiDownBBLinearStretched::~PhiPhiDownBBLinearStretched() {
    }

    void PhiPhiDownBBLinearStretched::operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim) {
      rec(source, result, index, dim, 0.0, 0.0);

    }

    void PhiPhiDownBBLinearStretched::rec(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim, float_t fl, float_t fr) {

      size_t seq = index.seq();

      float_t alpha_value = source[seq];

      SGPP::base::GridStorage::index_type::level_type current_level;
      SGPP::base::GridStorage::index_type::index_type current_index;

      index.get(dim, current_level, current_index);

      float_t posl = 0, posr = 0, currentPosition = 0;
      this->stretching->getAdjacentPositions(static_cast<int>(current_level), static_cast<int>(current_index), dim, currentPosition, posl, posr );
      float_t baseLength = posr - posl;
      float_t leftLength = currentPosition - posl;
      float_t rightLength = posr - currentPosition;

      // integration
      result[seq] = (1.0 / 3.0) * (baseLength) * alpha_value + fl / 6.0 * (baseLength + rightLength) + fr / 6.0 * (baseLength + leftLength);

      // dehierarchisation
      float_t fm  = (fr - fl) * (leftLength) / (baseLength) + fl + alpha_value;


      if (!index.hint()) {
        index.left_child(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, fl, fm);
        }

        index.step_right(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, fm, fr);
        }

        index.up(dim);
      }
    }

    // namespace detail

  } // namespace SGPP
}