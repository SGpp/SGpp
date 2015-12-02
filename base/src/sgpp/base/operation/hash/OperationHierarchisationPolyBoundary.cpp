/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)
#include <sgpp/base/operation/hash/OperationHierarchisationPolyBoundary.hpp>

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationPolyBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationPolyBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

namespace SGPP {
  namespace base {

    void OperationHierarchisationPolyBoundary::doHierarchisation(
      DataVector& node_values) {

      HierarchisationPolyBoundary func(this->storage, &this->base);
      sweep<HierarchisationPolyBoundary> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D_Boundary(node_values, node_values, i);
      }
    }

    void OperationHierarchisationPolyBoundary::doDehierarchisation(
      DataVector& alpha) {
      DehierarchisationPolyBoundary func(this->storage, &this->base);
      sweep<DehierarchisationPolyBoundary> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        DataVector source(alpha);
        s.sweep1D_Boundary(source, alpha, i);
      }
    }

  }
}
