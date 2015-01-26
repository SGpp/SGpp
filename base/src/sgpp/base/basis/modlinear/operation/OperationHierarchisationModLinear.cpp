/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)


#include "base/basis/modpoly/ModifiedPolyBasis.hpp"
#include "base/basis/modlinear/operation/OperationHierarchisationModLinear.hpp"
#include "base/basis/modlinear/algorithm_sweep/HierarchisationModLinear.hpp"
#include "base/basis/modlinear/algorithm_sweep/DehierarchisationModLinear.hpp"

#include "base/algorithm/sweep.hpp"


namespace sg {
  namespace base {
    /**
     * Implements the hierarchisation on a sprase grid with mod linear base functions
     *
     * @param node_values the functions values in the node base
     *
     * @todo (heinecke, nice) Implement the hierarchisation on the sparse grid with mod linear base functions
     */
    void OperationHierarchisationModLinear::doHierarchisation(DataVector& node_values) {
      HierarchisationModLinear func(this->storage);
      sweep<HierarchisationModLinear> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(node_values, node_values, i);
      }
    }

    /**
     * Implements the dehierarchisation on a sprase grid with mod linear base functions
     *
     * @param alpha the coefficients of the sparse grid's base functions
     *
     * @todo (heinecke, nice) Implement the dehierarchisation on the sparse grid with mod linear base functions
     */
    void OperationHierarchisationModLinear::doDehierarchisation(DataVector& alpha) {
      DehierarchisationModLinear func(this->storage);
      sweep<DehierarchisationModLinear> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(alpha, alpha, i);
      }
    }

  }
}
