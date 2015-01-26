/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/basis/linear/noboundary/operation/OperationHierarchisationLinear.hpp>
#include <sgpp/base/basis/linear/noboundary/algorithm_sweep/HierarchisationLinear.hpp>
#include <sgpp/base/basis/linear/noboundary/algorithm_sweep/DehierarchisationLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


namespace sg {
  namespace base {

    void OperationHierarchisationLinear::doHierarchisation(DataVector& node_values) {
      HierarchisationLinear func(this->storage);
      sweep<HierarchisationLinear> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(node_values, node_values, i);
      }
    }

    void OperationHierarchisationLinear::doDehierarchisation(DataVector& alpha) {
      DehierarchisationLinear func(this->storage);
      sweep<DehierarchisationLinear> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(alpha, alpha, i);
      }
    }

  }
}
