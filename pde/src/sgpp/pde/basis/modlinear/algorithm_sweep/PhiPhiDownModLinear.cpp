// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/modlinear/algorithm_sweep/PhiPhiDownModLinear.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {



    PhiPhiDownModLinear::PhiPhiDownModLinear(SGPP::base::GridStorage* storage) : storage(storage) {
    }

    PhiPhiDownModLinear::~PhiPhiDownModLinear() {
    }

    void PhiPhiDownModLinear::operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim) {
      rec(source, result, index, dim, 0.0, 0.0);
    }

    void PhiPhiDownModLinear::rec(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim, float_t fl, float_t fr) {
      size_t seq = index.seq();

      float_t alpha_value = source[seq];

      SGPP::base::GridStorage::index_type::level_type l;
      SGPP::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      float_t h = 1 / pow(2.0, static_cast<int>(l));
      float_t fm;

      // level 1, constant function
      if (l == 1) {
        //integration
        result[seq] = 0.0 + alpha_value;

        //dehierarchisation
        fm = (fl + fr) / 2.0 + alpha_value;

        //boundary value
        fl += alpha_value;
        fr += alpha_value;
      }
      // left boundary
      else if (i == 1) {
        //integration
        result[seq] = 2.0 / 3.0 * h * (2.0 * fl + fr)
                      + 8.0 / 3.0 * h * alpha_value;

        //dehierarchisation
        fm = (fl + fr) / 2.0 + alpha_value;

        //boundary value
        fl += 2.0 * alpha_value;
      }
      // right boundary
      else if (static_cast<int>(i) == static_cast<int>((1 << l) - 1)) {
        //integration
        result[seq] = 2.0 / 3.0 * h * (fl + 2.0 * fr)
                      + 8.0 / 3.0 * h * alpha_value;

        //dehierarchisation
        fm = (fl + fr) / 2.0 + alpha_value;

        //boundary value
        fr += 2.0 * alpha_value;
      }
      // inner functions
      else {
        //integration
        result[seq] = h * (fl + fr) / 2.0
                      + 2.0 / 3.0 * h * alpha_value;

        //dehierarchisation
        fm = (fl + fr) / 2.0 + alpha_value;

        //boundary value

      }

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



  }
}