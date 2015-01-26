/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Stefanie Schraufstetter (schraufs@in.tum.de)

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp>

namespace sg {
  namespace finance {



    DPhiPhiUpBBLinear::DPhiPhiUpBBLinear(sg::base::GridStorage* storage) : storage(storage), boundingBox(storage->getBoundingBox()) {
    }


    DPhiPhiUpBBLinear::~DPhiPhiUpBBLinear() {
    }

    void DPhiPhiUpBBLinear::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
      // get boundary values
      double fl = 0.0;
      double fr = 0.0;

      rec(source, result, index, dim, fl, fr);
    }


    void DPhiPhiUpBBLinear::rec(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr) {
      size_t seq = index.seq();

      fl = fr = 0.0;
      double fml = 0.0;
      double fmr = 0.0;

      sg::base::GridStorage::index_type::level_type current_level;
      sg::base::GridStorage::index_type::index_type current_index;

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

      double fm = fml + fmr;

      double alpha_value = source[seq];

      // transposed operations:
      result[seq] = fm;

      fl = (fm / 2.0) + (alpha_value * 0.5 + fl);
      fr = (fm / 2.0) + (alpha_value * (-0.5) + fr);
    }

    // namespace detail

  } // namespace sg
}
