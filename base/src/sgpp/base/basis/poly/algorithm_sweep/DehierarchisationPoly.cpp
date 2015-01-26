/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Kilian Roehner (roehner@tum.de)

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <cmath>

#include <sgpp/base/basis/poly/algorithm_sweep/DehierarchisationPoly.hpp>

namespace sg {

  namespace base {

    DehierarchisationPoly::DehierarchisationPoly(GridStorage* storage, SPolyBase* base) : storage(storage), base(base) {
    }

    DehierarchisationPoly::~DehierarchisationPoly() {
    }

    void DehierarchisationPoly::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim) {
      DataVector koeffs(index.getGridDepth(dim) + 1);
      koeffs.setAll(0.0);
      rec(source, result, index, dim, koeffs);
    }

    void DehierarchisationPoly::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, DataVector& koeffs) {
      // current position on the grid
      size_t seq = index.seq();

      level_type cur_lev;
      index_type cur_ind;

      // get current level and index from grid
      index.get(dim, cur_lev, cur_ind);

      // Dehierarchisation
      result[seq] = source[seq] + this->base->evalHierToTop(cur_lev, cur_ind, koeffs, cur_ind / (pow(2.0, static_cast<int>(cur_lev))));

      // recursive calls for the right and left side of the current node
      if (index.hint() == false) {
        koeffs[cur_lev] = result[seq];

        // descend left
        index.left_child(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, koeffs);
        }

        // descend right
        index.step_right(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, koeffs);
        }

        // ascend
        index.up(dim);

        koeffs[cur_lev] = 0.0;
      }
    }

  } // namespace base

} // namespace sg
