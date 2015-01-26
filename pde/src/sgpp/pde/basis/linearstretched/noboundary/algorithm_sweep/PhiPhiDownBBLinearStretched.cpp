/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp>


namespace sg {
  namespace pde {



    PhiPhiDownBBLinearStretched::PhiPhiDownBBLinearStretched(sg::base::GridStorage* storage) : storage(storage), stretching(storage->getStretching()) {
    }

    PhiPhiDownBBLinearStretched::~PhiPhiDownBBLinearStretched() {
    }

    void PhiPhiDownBBLinearStretched::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
      rec(source, result, index, dim, 0.0, 0.0);

    }

    void PhiPhiDownBBLinearStretched::rec(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double fl, double fr) {

      size_t seq = index.seq();

      double alpha_value = source[seq];

      sg::base::GridStorage::index_type::level_type current_level;
      sg::base::GridStorage::index_type::index_type current_index;

      index.get(dim, current_level, current_index);

      double posl = 0, posr = 0, currentPosition = 0;
      this->stretching->getAdjacentPositions(static_cast<int>(current_level), static_cast<int>(current_index), dim, currentPosition, posl, posr );
      double baseLength = posr - posl;
      double leftLength = currentPosition - posl;
      double rightLength = posr - currentPosition;

      // integration
      result[seq] = (1.0 / 3.0) * (baseLength) * alpha_value + fl / 6.0 * (baseLength + rightLength) + fr / 6.0 * (baseLength + leftLength);

      // dehierarchisation
      double fm  = (fr - fl) * (leftLength) / (baseLength) + fl + alpha_value;


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

  } // namespace sg
}
