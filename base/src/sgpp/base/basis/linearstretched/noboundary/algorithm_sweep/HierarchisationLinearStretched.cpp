/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/basis/linearstretched/noboundary/algorithm_sweep/HierarchisationLinearStretched.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {



    HierarchisationLinearStretched::HierarchisationLinearStretched(GridStorage* storage) : storage(storage), stretch(storage->getStretching()) {
    }

    HierarchisationLinearStretched::~HierarchisationLinearStretched() {
    }

    void HierarchisationLinearStretched::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim) {
      rec(source, result, index, dim, 0.0, 0.0);
    }

    void HierarchisationLinearStretched::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr) {
      // current position on the grid
      size_t seq = index.seq();
      // value in the middle, needed for recursive call and calculation of the hierarchical surplus
      double fm = source[seq];

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

      // hierarchisation
      GridStorage::index_type::level_type current_level;
      GridStorage::index_type::index_type current_index;
      index.get(dim, current_level, current_index);

      double posl = 0, posr = 0, posc = 0;

      if ((static_cast<int>(current_level)) == 0) {
        std::cout << "printing fl and fr " << fl << " " << fr << std::endl;
      }


      stretch->getAdjacentPositions(static_cast<int>(current_level), static_cast<int>(current_index), dim, posc, posl, posr );


      double fcurr = (fr - fl) * (posc - posl) / (posr - posl) + fl;
      result[seq] = fm - fcurr;

    }

    // namespace detail

  } // namespace SGPP
}
