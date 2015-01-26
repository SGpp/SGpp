/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)



#include <sgpp/base/basis/linearstretched/noboundary/operation/OperationHierarchisationLinearStretched.hpp>
#include <sgpp/base/basis/linearstretched/noboundary/algorithm_sweep/HierarchisationLinearStretched.hpp>
#include <sgpp/base/basis/linearstretched/noboundary/algorithm_sweep/DehierarchisationLinearStretched.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


namespace sg {
  namespace base {

    void OperationHierarchisationLinearStretched::doHierarchisation(DataVector& node_values) {
      HierarchisationLinearStretched func(this->storage);
      sweep<HierarchisationLinearStretched> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(node_values, node_values, i);
      }
    }

    void OperationHierarchisationLinearStretched::doDehierarchisation(DataVector& alpha) {
      DehierarchisationLinearStretched func(this->storage);
      sweep<DehierarchisationLinearStretched> s(func, this->storage);

      // Execute hierarchisation in every dimension of the grid
      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(alpha, alpha, i);
      }
    }

  }
}
