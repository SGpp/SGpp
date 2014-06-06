/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "finance/basis/linear/noboundary/algorithm_sweep/XPhidPhiDownBBLinear.hpp"

namespace sg {
  namespace finance {



    XPhidPhiDownBBLinear::XPhidPhiDownBBLinear(sg::base::GridStorage* storage) : storage(storage), boundingBox(storage->getBoundingBox()) {
    }

    XPhidPhiDownBBLinear::~XPhidPhiDownBBLinear() {
    }

    void XPhidPhiDownBBLinear::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
      double q = boundingBox->getIntervalWidth(dim);
      double t = boundingBox->getIntervalOffset(dim);

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

    void XPhidPhiDownBBLinear::rec(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double fl, double fr) {
      size_t seq = index.seq();

      double alpha_value = source[seq];

      sg::base::GridStorage::index_type::level_type l;
      sg::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      double hhalf = 1.0 / static_cast<double>(1 << (l + 1));
      double i_dbl = static_cast<double>(i);

      // integration
      result[seq] = (  ( (fl * ((hhalf * i_dbl) - hhalf)) + (fr * (((-1.0) * (hhalf * i_dbl)) - hhalf)) ) - ((1.0 / 3.0) * (((1.0 / static_cast<double>(1 << l))) * alpha_value)) ); // diagonal entry

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

    void XPhidPhiDownBBLinear::recBB(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double fl, double fr, double q, double t) {
      size_t seq = index.seq();

      double alpha_value = source[seq];

      sg::base::GridStorage::index_type::level_type l;
      sg::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      double hhalf = 1.0 / static_cast<double>(1 << (l + 1));
      double i_dbl = static_cast<double>(i);

      // integration
      result[seq] = (  ( (fl * ((q * ((hhalf * i_dbl) - hhalf)) + (0.5 * t))) + (fr * ((q * (((-1.0) * (hhalf * i_dbl)) - hhalf)) - (0.5 * t))) )
                       - ((1.0 / 3.0) * (((1.0 / static_cast<double>(1 << l)) * q) * alpha_value)) ); // diagonal entry

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
