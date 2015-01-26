// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/basis/linearstretched/boundary/operation/OperationHierarchisationLinearStretchedBoundary.hpp>
#include <sgpp/base/basis/linearstretched/boundary/algorithm_sweep/HierarchisationLinearStretchedBoundary.hpp>
#include <sgpp/base/basis/linearstretched/boundary/algorithm_sweep/DehierarchisationLinearStretchedBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    void OperationHierarchisationLinearStretchedBoundary::doHierarchisation(DataVector& node_values) {
      HierarchisationLinearStretchedBoundary func(this->storage);
      sweep<HierarchisationLinearStretchedBoundary> s(func, this->storage);

      // N D case
      if (this->storage->dim() > 1) {
        for (size_t i = 0; i < this->storage->dim(); i++) {
          s.sweep1D_Boundary(node_values, node_values, i);
        }
      }
      // 1 D case
      else {
        s.sweep1D(node_values, node_values, 0);
      }
    }

    void OperationHierarchisationLinearStretchedBoundary::doDehierarchisation(DataVector& alpha) {
      DehierarchisationLinearStretchedBoundary func(this->storage);
      sweep<DehierarchisationLinearStretchedBoundary> s(func, this->storage);

      // N D case
      if (this->storage->dim() > 1) {
        for (size_t i = 0; i < this->storage->dim(); i++) {
          s.sweep1D_Boundary(alpha, alpha, i);
        }
      }
      // 1 D case
      else {
        s.sweep1D(alpha, alpha, 0);
      }
    }

  }
}