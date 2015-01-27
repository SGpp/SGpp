// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/StencilDehierarchisationLinear.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {



    StencilDehierarchisationLinear::StencilDehierarchisationLinear(
      GridStorage* storage,
      OperationStencilHierarchisation::IndexStencil& surplusStencil,
      OperationStencilHierarchisation::IndexStencil& neighborStencil,
      OperationStencilHierarchisation::WeightStencil& weightStencil) :
      storage(storage), _surplusStencil(surplusStencil),
      _neighborStencil(neighborStencil), _weightStencil(weightStencil) {
    }

    StencilDehierarchisationLinear::~StencilDehierarchisationLinear() {
    }

    void StencilDehierarchisationLinear::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim) {
      rec(source, result, index, dim, -1, -1);
    }

    void StencilDehierarchisationLinear::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, int seql, int seqr) {
      // current position on the grid
      int seqm = (int) index.seq();

      // dehierarchisation
      if (seql >= 0) {
        _surplusStencil.push_back(seqm);
        _neighborStencil.push_back(seql);
        _weightStencil.push_back(0.5f);
      }

      if (seqr >= 0) {
        _surplusStencil.push_back(seqm);
        _neighborStencil.push_back(seqr);
        _weightStencil.push_back(0.5f);
      }

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
    }

    // namespace detail

  } // namespace SGPP
}