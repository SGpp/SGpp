// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationLinear.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationLinear.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
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