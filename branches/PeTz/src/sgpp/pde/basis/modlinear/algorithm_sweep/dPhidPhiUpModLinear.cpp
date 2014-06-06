/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/basis/modlinear/algorithm_sweep/dPhidPhiUpModLinear.hpp"

namespace sg {
  namespace pde {



    dPhidPhiUpModLinear::dPhidPhiUpModLinear(sg::base::GridStorage* storage) : storage(storage) {
    }

    dPhidPhiUpModLinear::~dPhidPhiUpModLinear() {
    }

    void dPhidPhiUpModLinear::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
      double f = 0.0;
      rec(source, result, index, dim, f);
    }

    void dPhidPhiUpModLinear::rec(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double& f) {
      size_t seq = index.seq();

      sg::base::GridStorage::index_type::level_type l;
      sg::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      double alpha_value = source[seq];
      double ht = pow(2.0, static_cast<int>(l));

      if (l == 1) {
        f = 0.0;

        if (!index.hint()) {
          index.left_child(dim);

          if (!storage->end(index.seq())) {
            rec(source, result, index, dim, f);
          }

          f = 0.0;
          index.step_right(dim);

          if (!storage->end(index.seq())) {
            rec(source, result, index, dim, f);
          }

          index.up(dim);
        }

        result[seq] = 0.0;
      }
      // left boundary
      else if (i == 1) {
        f = 0.0;

        if (!index.hint()) {
          index.left_child(dim);

          if (!storage->end(index.seq())) {
            rec(source, result, index, dim, f);
          }

          index.up(dim);
        }

        result[seq] = ht * f;

        f += 2.0 * alpha_value;
      }
      // right boundary
      else if (static_cast<int>(i) == static_cast<int>((1 << l) - 1)) {
        f = 0.0;

        if (!index.hint()) {
          index.right_child(dim);

          if (!storage->end(index.seq())) {
            rec(source, result, index, dim, f);
          }

          index.up(dim);
        }

        result[seq] = ht * f;

        f += 2.0 * alpha_value;
      }
    }



  }
}
