// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/XdPhiPhiDownBBLinearStretched.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {



    XdPhiPhiDownBBLinearStretched::XdPhiPhiDownBBLinearStretched(SGPP::base::GridStorage* storage) : storage(storage), stretching(storage->getStretching()) {
    }

    XdPhiPhiDownBBLinearStretched::~XdPhiPhiDownBBLinearStretched() {
    }

    void XdPhiPhiDownBBLinearStretched::operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim) {

      rec(source, result, index, dim, 0.0, 0.0);

    }


    void XdPhiPhiDownBBLinearStretched::rec(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim, float_t fl, float_t fr) {
      size_t seq = index.seq();

      float_t alpha_value = source[seq];

      SGPP::base::GridStorage::index_type::level_type l;
      SGPP::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);
      //get the positions of the current index as well as its left and right neighbors
      float_t posl = 0, posr = 0, posc = 0;
      this->stretching->getAdjacentPositions(static_cast<int>(l), static_cast<int>(i), dim, posc, posl, posr );

      float_t baseLength = posr - posl;
      float_t leftLength = posc - posl;

      result[seq] = fl * (-1.0 / 6.0) * (posc + posl + posr) - fr * (-1.0 / 6.0) * (posc + posl + posr)
                    - (1.0 / 6.0) * baseLength * alpha_value; // diagonal entry




      // dehierarchisation
      float_t fm =  (fr - fl) * (leftLength) / (baseLength) + fl + alpha_value;


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

  } // namespace SGPP
}