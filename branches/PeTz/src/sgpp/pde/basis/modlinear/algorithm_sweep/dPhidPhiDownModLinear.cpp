/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/basis/modlinear/algorithm_sweep/dPhidPhiDownModLinear.hpp"

namespace sg {
  namespace pde {



    dPhidPhiDownModLinear::dPhidPhiDownModLinear(sg::base::GridStorage* storage) : storage(storage) {
    }

    dPhidPhiDownModLinear::~dPhidPhiDownModLinear() {
    }

    void dPhidPhiDownModLinear::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
      rec(source, result, index, dim, 0.0);
    }

    void dPhidPhiDownModLinear::rec(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double f) {
      size_t seq = index.seq();
      sg::base::GridStorage::index_type::level_type l;
      sg::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      double alpha_value = source[seq];
      double ht = pow(2.0, static_cast<int>(l));
      double f_local = 0.0;

      // level 1, constant function
      if (l == 1) {
        f_local = 0.0;
        result[seq] = 0.0
                      + 0.0;
      }
      // left boundary & right boundary
      else if ((i == 1) || (static_cast<int>(i) == static_cast<int>((1 << l) - 1))) {
        f_local = ht * alpha_value;
        result[seq] = 2.0 * f
                      + 2.0 * f_local;
      }
      // inner functions
      else {
        f_local = ht * alpha_value;
        result[seq] = 0.0
                      + 2.0 * f_local;
      }

      if (!index.hint()) {
        index.left_child(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, f + f_local);
        }

        index.step_right(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, f + f_local);
        }

        index.up(dim);
      }

    }



  }
}
