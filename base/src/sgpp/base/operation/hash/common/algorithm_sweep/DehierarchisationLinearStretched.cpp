// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationLinearStretched.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {



    DehierarchisationLinearStretched::DehierarchisationLinearStretched(GridStorage* storage) : storage(storage) {
    }

    DehierarchisationLinearStretched::~DehierarchisationLinearStretched() {
    }

    void DehierarchisationLinearStretched::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim) {
      rec(source, result, index, dim, 0.0, 0.0);
    }

    void DehierarchisationLinearStretched::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr) {
      // current position on the grid
      size_t seq = index.seq();
      // value in the middle, needed for recursive call and calculation of the hierarchical surplus
      double fm = source[seq];

      // dehierarchisation
      GridStorage::index_type::level_type current_level;
      GridStorage::index_type::index_type current_index;
      index.get(dim, current_level, current_index);

      //get the positions of the current index as well as its left and right neighbors
      double posl = 0, posr = 0, posc = 0;
      storage->getStretching()->getAdjacentPositions(static_cast<int>(current_level), static_cast<int>(current_index), dim, posc, posl, posr );

      double fcurr = (fr - fl) * (posc - posl) / (posr - posl) + fl;
      fm += fcurr;
      result[seq] = fm;


      // recursive calls for the right and left side of the current node
      if (index.hint() == false) {
        // descend left
        index.left_child(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, fl, fm);
        }

        // descend right
        index.step_right(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, fm, fr);
        }

        // ascend
        index.up(dim);
      }
    }

    // namespace detail

  } // namespace SGPP
}