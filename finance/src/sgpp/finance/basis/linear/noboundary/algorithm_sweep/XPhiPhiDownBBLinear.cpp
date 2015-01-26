/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de) Chao qi (qic@in.tum.de) Stefanie Schraufstetter (schraufs@in.tum.de)

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp>

namespace sg {
  namespace finance {



    XPhiPhiDownBBLinear::XPhiPhiDownBBLinear(sg::base::GridStorage* storage) : storage(storage), boundingBox(storage->getBoundingBox()) {
    }

    XPhiPhiDownBBLinear::~XPhiPhiDownBBLinear() {
    }

    void XPhiPhiDownBBLinear::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
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

    void XPhiPhiDownBBLinear::rec(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double fl, double fr) {
      size_t seq = index.seq();

      double alpha_value = source[seq];

      sg::base::GridStorage::index_type::level_type l;
      sg::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      double i_dbl = static_cast<double>(i);
      int l_int = static_cast<int>(l);

      double hsquare = (1.0 / static_cast<double>(1 << (2 * l_int)));

      // integration
      result[seq] = (hsquare * ((fl + fr) / 2.0)) * i_dbl + hsquare * (fr - fl) / 12.0 + (((2.0 / 3.0) * hsquare * i_dbl) * alpha_value);

      // dehierarchisation
      double fm = ((fl + fr) / 2.0) + alpha_value;

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

    void XPhiPhiDownBBLinear::recBB(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double fl, double fr, double q, double t) {
      size_t seq = index.seq();

      double alpha_value = source[seq];

      sg::base::GridStorage::index_type::level_type l;
      sg::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      double i_dbl = static_cast<double>(i);
      int l_int = static_cast<int>(l);

      double h = (1.0 / (static_cast<double>(1 << (l_int))));

      // integration
      result[seq] = (h * h * i_dbl * q * q + h * q * t) * ((fl + fr) / 2.0) + h * h * q * q * (fr - fl) / 12.0 + (((2.0 / 3.0) * h * q * t + (2.0 / 3.0) * i_dbl * h * h * q * q ) * alpha_value); // diagonal entry

      // dehierarchisation
      double fm = ((fl + fr) / 2.0) + alpha_value;

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

  } // namespace sg
}
