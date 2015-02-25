// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/PhidPhiDownBBLinearStretched.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {



    PhidPhiDownBBLinearStretched::PhidPhiDownBBLinearStretched(SGPP::base::GridStorage* storage) : storage(storage), stretching(storage->getStretching()) {
    }

    PhidPhiDownBBLinearStretched::~PhidPhiDownBBLinearStretched() {
    }

    void PhidPhiDownBBLinearStretched::operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim) {
      rec(source, result, index, dim, 0.0, 0.0);
    }

    void PhidPhiDownBBLinearStretched::rec(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim, float_t fl, float_t fr) {
      //Fix the eqn (Selcuk)
      size_t seq = index.seq();

      float_t alpha_value = source[seq];

      SGPP::base::GridStorage::index_type::level_type l;
      SGPP::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      float_t posl = 0, posr = 0, posc = 0;
      this->stretching->getAdjacentPositions(static_cast<int>(l), static_cast<int>(i), dim, posc, posl, posr );

      // integration
      result[seq] = (  0.5 * (fl - fr) ); // diagonal entry = 0.0

      // dehierarchisation
      float_t fm  = (fr - fl) * (posc - posl) / (posr - posl) + fl + alpha_value;

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