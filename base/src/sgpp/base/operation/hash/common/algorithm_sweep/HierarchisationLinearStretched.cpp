// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationLinearStretched.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {



HierarchisationLinearStretched::HierarchisationLinearStretched(
  GridStorage* storage) : storage(storage), stretch(storage->getStretching()) {
}

HierarchisationLinearStretched::~HierarchisationLinearStretched() {
}

void HierarchisationLinearStretched::operator()(DataVector& source,
    DataVector& result, grid_iterator& index, size_t dim) {
  rec(source, result, index, dim, 0.0, 0.0);
}

void HierarchisationLinearStretched::rec(DataVector& source, DataVector& result,
    grid_iterator& index, size_t dim, float_t fl, float_t fr) {
  // current position on the grid
  size_t seq = index.seq();
  // value in the middle, needed for recursive call and calculation of the hierarchical surplus
  float_t fm = source[seq];

  // recursive calls for the right and left side of the current node
  if (index.hint() == false) {
    // descend left
    index.leftChild(dim);

    if (!storage->end(index.seq())) {
      rec(source, result, index, dim, fl, fm);
    }

    // descend right
    index.stepRight(dim);

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

  float_t posl = 0, posr = 0, posc = 0;

  if ((static_cast<int>(current_level)) == 0) {
    std::cout << "printing fl and fr " << fl << " " << fr << std::endl;
  }


  stretch->getAdjacentPositions(static_cast<int>(current_level),
                                static_cast<int>(current_index), dim, posc, posl, posr );


  float_t fcurr = (fr - fl) * (posc - posl) / (posr - posl) + fl;
  result[seq] = fm - fcurr;

}

// namespace detail

} // namespace SGPP
}