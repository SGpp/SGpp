/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/basis/modlinear/algorithm_sweep/PhiPhiDownModLinear.hpp"

namespace sg {
  namespace pde {



    PhiPhiDownModLinear::PhiPhiDownModLinear(sg::base::GridStorage* storage) : storage(storage) {
    }

    PhiPhiDownModLinear::~PhiPhiDownModLinear() {
    }

    void PhiPhiDownModLinear::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {
      rec(source, result, index, dim, 0.0, 0.0);
    }

    void PhiPhiDownModLinear::rec(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double fl, double fr) {
      size_t seq = index.seq();

      double alpha_value = source[seq];

      sg::base::GridStorage::index_type::level_type l;
      sg::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      double h = 1 / pow(2.0, static_cast<int>(l));
      double fm;

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
