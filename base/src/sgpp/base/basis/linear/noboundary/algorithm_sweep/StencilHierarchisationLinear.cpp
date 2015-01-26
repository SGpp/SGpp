/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Gerrit Buse (buse@in.tum.de)

#include <iostream>

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/basis/linear/noboundary/algorithm_sweep/StencilHierarchisationLinear.hpp>

namespace sg {
  namespace base {



    StencilHierarchisationLinear::StencilHierarchisationLinear(
      GridStorage* storage,
      OperationStencilHierarchisation::IndexStencil& surplusStencil,
      OperationStencilHierarchisation::IndexStencil& neighborStencil,
      OperationStencilHierarchisation::WeightStencil& weightStencil) :
      storage(storage), _surplusStencil(surplusStencil),
      _neighborStencil(neighborStencil), _weightStencil(weightStencil) {
    }

    StencilHierarchisationLinear::~StencilHierarchisationLinear() {
    }

    void StencilHierarchisationLinear::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim) {
      rec(source, result, index, dim, -1, -1);
    }

    void StencilHierarchisationLinear::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, int seql, int seqr) {
      // current position on the grid
      int seqm = (int) index.seq();

      // recursive calls for the right and left side of the current node
      if (index.hint() == false) {
        // descend left
        index.left_child(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, seql, seqm);
        }

        // descend right
        index.step_right(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, seqm, seqr);
        }

        // ascend
        index.up(dim);
      }

      // hierarchisation

      if (seql >= 0) {
        _surplusStencil.push_back(seqm);
        _neighborStencil.push_back(seql);
        _weightStencil.push_back(-0.5f);
      }

      if (seqr >= 0) {
        _surplusStencil.push_back(seqm);
        _neighborStencil.push_back(seqr);
        _weightStencil.push_back(-0.5f);
      }
    }

    // namespace detail

  } // namespace sg
}

