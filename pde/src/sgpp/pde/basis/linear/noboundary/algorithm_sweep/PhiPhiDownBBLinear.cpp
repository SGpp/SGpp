/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"

namespace sg {
  namespace pde {



    PhiPhiDownBBLinear::PhiPhiDownBBLinear(sg::base::GridStorage* storage) : storage(storage), boundingBox(storage->getBoundingBox()) {
    }

    PhiPhiDownBBLinear::~PhiPhiDownBBLinear() {
    }

    void PhiPhiDownBBLinear::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
      //std::cout << dim << std::endl;
      //std::cout << index.toString() << std::endl;

      double q = this->boundingBox->getIntervalWidth(dim);
      double t = this->boundingBox->getIntervalOffset(dim);

      bool useBB = false;

      if (q != 1.0 || t != 0.0) {
        useBB = true;
      }

      if (useBB) {
        recBB(source, result, index, dim, 0.0, 0.0, q, t);
      } else {
        rec(source, result, index, dim, 0.0, 0.0);
      }
    }

    void PhiPhiDownBBLinear::rec(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double fl, double fr) {
      size_t seq = index.seq();

      double alpha_value = source[seq];

      sg::base::GridStorage::index_type::level_type l;
      sg::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      double h = 1.0 / static_cast<double>(1 << l);
      double tmp_m = ((fl + fr) / 2.0);

      // integration
      result[seq] = (h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value);

      // dehierarchisation
      double fm = tmp_m + alpha_value;

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

    void PhiPhiDownBBLinear::recBB(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double fl, double fr, double q, double t) {
      size_t seq = index.seq();

      double alpha_value = source[seq];

      sg::base::GridStorage::index_type::level_type l;
      sg::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      double h = 1.0 / static_cast<double>(1 << l);
      double tmp_m = ((fl + fr) / 2.0);

      // integration
      result[seq] = ((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value)) * q;  // diagonal entry

      // dehierarchisation
      double fm = tmp_m + alpha_value;

      if (!index.hint()) {
        index.left_child(dim);

        if (!storage->end(index.seq())) {
          recBB(source, result, index, dim, fl, fm, q, t);
        }

        index.step_right(dim);

        if (!storage->end(index.seq())) {
          recBB(source, result, index, dim, fm, fr, q, t);
        }

        index.up(dim);
      }
    }

    // namespace detail
  }
  // namespace sg
}
