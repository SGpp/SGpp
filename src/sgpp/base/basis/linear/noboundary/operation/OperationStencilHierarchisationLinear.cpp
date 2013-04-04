/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Gerrit Buse (buse@in.tum.de)

#include "base/basis/linear/noboundary/operation/OperationStencilHierarchisationLinear.hpp"
#include "base/basis/linear/noboundary/algorithm_sweep/StencilHierarchisationLinear.hpp"
#include "base/basis/linear/noboundary/algorithm_sweep/StencilDehierarchisationLinear.hpp"

#include "base/algorithm/sweep.hpp"


namespace sg {
  namespace base {

    void OperationStencilHierarchisationLinear::doHierarchisation(DataVector& node_values) {
      surplusStencil.clear();
      neighborStencil.clear();
      weightStencil.clear();
      StencilHierarchisationLinear func(this->storage, surplusStencil, neighborStencil, weightStencil);
      sweep<StencilHierarchisationLinear> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(node_values, node_values, i);
      }
    }

    void OperationStencilHierarchisationLinear::doDehierarchisation(DataVector& alpha) {
      surplusStencil.clear();
      neighborStencil.clear();
      weightStencil.clear();
      StencilDehierarchisationLinear func(this->storage, surplusStencil, neighborStencil, weightStencil);
      sweep<StencilDehierarchisationLinear> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(alpha, alpha, i);
      }
    }

  }
}
