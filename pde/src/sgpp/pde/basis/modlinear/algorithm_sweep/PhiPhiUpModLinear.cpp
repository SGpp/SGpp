// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/modlinear/algorithm_sweep/PhiPhiUpModLinear.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {



    PhiPhiUpModLinear::PhiPhiUpModLinear(SGPP::base::GridStorage* storage) : storage(storage) {
    }

    PhiPhiUpModLinear::~PhiPhiUpModLinear() {
    }

    void PhiPhiUpModLinear::operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim) {
      float_t fl = 0.0;
      float_t fr = 0.0;
      rec(source, result, index, dim, fl, fr);
    }

    void PhiPhiUpModLinear::rec(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim, float_t& fl, float_t& fr) {
      size_t seq = index.seq();

      float_t alpha_value = source[seq];

      SGPP::base::GridStorage::index_type::level_type l;
      SGPP::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      float_t h = 1 / pow(2.0, static_cast<int>(l));
      float_t fml = 0.0;
      float_t fmr = 0.0;

      if (!index.hint()) {
        index.leftChild(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, fl, fml);
        }

        index.stepRight(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, fmr, fr);
        }

        index.up(dim);
      }

      float_t fm = fml + fmr;

      // level 1, constant function
      if (l == 1) {
        result[seq] = fl + fm + fr;

        fl += fm / 2.0 + alpha_value;
        fr += fm / 2.0 + alpha_value;
      }
      // left boundary
      else if (i == 1) {
        result[seq] = 2.0 * fl + fm;

        fl += fm / 2.0 + 4.0 / 3.0 * h * alpha_value;
        fr += fm / 2.0 + 2.0 / 3.0 * h * alpha_value;
      }
      // right boundary
      else if (static_cast<int>(i) == static_cast<int>((1 << l) - 1)) {
        result[seq] = 2.0 * fr + fm;

        fl += fm / 2.0 + 2.0 / 3.0 * h * alpha_value;
        fr += fm / 2.0 + 4.0 / 3.0 * h * alpha_value;
      }
      // inner functions
      else {
        result[seq] = fm;

        fl += fm / 2.0 + h / 2.0 * alpha_value;
        fr += fm / 2.0 + h / 2.0 * alpha_value;
      }
    }



  }
}