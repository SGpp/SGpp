/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "finance/basis/linearstretched/noboundary/algorithm_sweep/XdPhiPhiDownBBLinearStretched.hpp"

namespace sg {
  namespace finance {



    XdPhiPhiDownBBLinearStretched::XdPhiPhiDownBBLinearStretched(sg::base::GridStorage* storage) : storage(storage), stretching(storage->getStretching()) {
    }

    XdPhiPhiDownBBLinearStretched::~XdPhiPhiDownBBLinearStretched() {
    }

    void XdPhiPhiDownBBLinearStretched::operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim) {

      rec(source, result, index, dim, 0.0, 0.0);

    }


    void XdPhiPhiDownBBLinearStretched::rec(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double fl, double fr) {
      size_t seq = index.seq();

      double alpha_value = source[seq];

      sg::base::GridStorage::index_type::level_type l;
      sg::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);
      //get the positions of the current index as well as its left and right neighbors
      double posl = 0, posr = 0, posc = 0;
      this->stretching->getAdjacentPositions(static_cast<int>(l), static_cast<int>(i), dim, posc, posl, posr );

      double baseLength = posr - posl;
      double leftLength = posc - posl;

      result[seq] = fl * (-1.0 / 6.0) * (posc + posl + posr) - fr * (-1.0 / 6.0) * (posc + posl + posr)
                    - (1.0 / 6.0) * baseLength * alpha_value; // diagonal entry




      // dehierarchisation
      double fm =  (fr - fl) * (leftLength) / (baseLength) + fl + alpha_value;


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

    // namespace detail

  } // namespace sg
}
