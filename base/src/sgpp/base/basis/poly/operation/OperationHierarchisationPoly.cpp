/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/basis/poly/operation/OperationHierarchisationPoly.hpp>

#include <sgpp/base/basis/poly/algorithm_sweep/HierarchisationPoly.hpp>
#include <sgpp/base/basis/poly/algorithm_sweep/DehierarchisationPoly.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    void OperationHierarchisationPoly::doHierarchisation(DataVector& node_values) {

      HierarchisationPoly func(this->storage, &this->base);
      sweep<HierarchisationPoly> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(node_values, node_values, i);
      }
    }

    void OperationHierarchisationPoly::doDehierarchisation(DataVector& alpha) {
      DehierarchisationPoly func(this->storage, &this->base);
      sweep<DehierarchisationPoly> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(alpha, alpha, i);
      }
    }

  }
}
