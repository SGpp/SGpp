/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (dirk.pflueger@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

namespace sg {
  namespace pde {



    PhiPhiUpBBLinear::PhiPhiUpBBLinear(sg::base::GridStorage* storage) : storage(storage), boundingBox(storage->getBoundingBox()) {
    }

    PhiPhiUpBBLinear::~PhiPhiUpBBLinear() {
    }

    void PhiPhiUpBBLinear::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
      double q = boundingBox->getIntervalWidth(dim);
      double t = boundingBox->getIntervalOffset(dim);

      bool useBB = false;

      if (q != 1.0 || t != 0.0) {
        useBB = true;
      }

      // get boundary values
      double fl = 0.0;
      double fr = 0.0;

      if (useBB) {
        recBB(source, result, index, dim, fl, fr, q, t);
      } else {
        rec(source, result, index, dim, fl, fr);
      }
    }

    void PhiPhiUpBBLinear::rec(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr) {
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

      double tmp = (fm / 2.0) + (alpha_value / static_cast<double>(1 << (current_level + 1)));

      fl = tmp + fl;
      fr = tmp + fr;
    }

    void PhiPhiUpBBLinear::recBB(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr, double q, double t) {
      size_t seq = index.seq();

      fl = fr = 0.0;
      double fml = 0.0;
      double fmr = 0.0;

      sg::base::GridStorage::index_type::level_type current_level;
      sg::base::GridStorage::index_type::index_type current_index;

      if (!index.hint()) {
        index.left_child(dim);

        if (!storage->end(index.seq())) {
          recBB(source, result, index, dim, fl, fml, q, t);
        }

        index.step_right(dim);

        if (!storage->end(index.seq())) {
          recBB(source, result, index, dim, fmr, fr, q , t);
        }

        index.up(dim);
      }

      index.get(dim, current_level, current_index);

      double fm = fml + fmr;
      double alpha_value = source[seq];

      // transposed operations:
      result[seq] = fm;

      double tmp = ((fm / 2.0) + ((alpha_value / static_cast<double>(1 << (current_level + 1))) * q));

      fl = tmp + fl;
      fr = tmp + fr;
    }

    // namespace detail
  }
  // namespace sg
}
