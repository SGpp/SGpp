// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/SqXdPhidPhiUpBBLinearStretched.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {



    SqXdPhidPhiUpBBLinearStretched::SqXdPhidPhiUpBBLinearStretched(SGPP::base::GridStorage* storage) : storage(storage), stretching(storage->getStretching()) {
    }


    SqXdPhidPhiUpBBLinearStretched::~SqXdPhidPhiUpBBLinearStretched() {
    }


    void SqXdPhidPhiUpBBLinearStretched::operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim) {

      // get boundary values
      double fl = 0.0;
      double fr = 0.0;

      rec(source, result, index, dim, fl, fr);

    }

    void SqXdPhidPhiUpBBLinearStretched::rec(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr) {
      size_t seq = index.seq();

      fl = fr = 0.0;
      double fml = 0.0;
      double fmr = 0.0;

      SGPP::base::GridStorage::index_type::level_type current_level;
      SGPP::base::GridStorage::index_type::index_type current_index;

      if (!index.hint()) {
        index.left_child(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, fl, fml);
        }

        index.step_right(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, fmr, fr);
        }

        index.up(dim);
      }

      index.get(dim, current_level, current_index);
      double posl = 0, posr = 0, posc = 0;
      this->stretching->getAdjacentPositions(static_cast<int>(current_level), static_cast<int>(current_index), dim, posc, posl, posr);
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


    // namespace detail

  } // namespace SGPP
}