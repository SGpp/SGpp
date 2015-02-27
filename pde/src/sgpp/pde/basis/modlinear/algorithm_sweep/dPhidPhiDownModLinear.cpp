// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/modlinear/algorithm_sweep/dPhidPhiDownModLinear.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {



    dPhidPhiDownModLinear::dPhidPhiDownModLinear(SGPP::base::GridStorage* storage) : storage(storage) {
    }

    dPhidPhiDownModLinear::~dPhidPhiDownModLinear() {
    }

    void dPhidPhiDownModLinear::operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim) {
      rec(source, result, index, dim, 0.0);
    }

    void dPhidPhiDownModLinear::rec(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim, float_t f) {
      size_t seq = index.seq();
      SGPP::base::GridStorage::index_type::level_type l;
      SGPP::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      float_t alpha_value = source[seq];
      float_t ht = pow(2.0, static_cast<int>(l));
      float_t f_local = 0.0;

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