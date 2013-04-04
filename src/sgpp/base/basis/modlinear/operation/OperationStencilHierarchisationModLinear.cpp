/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Gerrit Buse (buse@in.tum.de)

#include "base/basis/modlinear/operation/OperationStencilHierarchisationModLinear.hpp"
#include "base/basis/modlinear/algorithm_sweep/StencilHierarchisationModLinear.hpp"
#include "base/basis/modlinear/algorithm_sweep/StencilDehierarchisationModLinear.hpp"

#include "base/algorithm/sweep.hpp"


namespace sg {
  namespace base {

    void OperationStencilHierarchisationModLinear::doHierarchisation(DataVector& node_values) {
      surplusStencil.clear();
      neighborStencil.clear();
      weightStencil.clear();
      StencilHierarchisationModLinear func(this->storage, surplusStencil, neighborStencil, weightStencil);
      sweep<StencilHierarchisationModLinear> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(node_values, node_values, i);
      }
    }

    void OperationStencilHierarchisationModLinear::doDehierarchisation(DataVector& alpha) {
      surplusStencil.clear();
      neighborStencil.clear();
      weightStencil.clear();
      StencilDehierarchisationModLinear func(this->storage, surplusStencil, neighborStencil, weightStencil);
      sweep<StencilDehierarchisationModLinear> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(alpha, alpha, i);
      }
    }

  }
}
