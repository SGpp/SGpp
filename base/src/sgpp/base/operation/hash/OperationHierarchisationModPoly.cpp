// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationModPoly.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationModPoly.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationModPoly.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    void OperationHierarchisationModPoly::doHierarchisation(DataVector& node_values) {
      HierarchisationModPoly func(this->storage, &this->base);
      sweep<HierarchisationModPoly> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(node_values, node_values, i);
      }
    }

    void OperationHierarchisationModPoly::doDehierarchisation(DataVector& alpha) {
      DehierarchisationModPoly func(this->storage, &this->base);
      sweep<DehierarchisationModPoly> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(alpha, alpha, i);
      }
    }

  }
}