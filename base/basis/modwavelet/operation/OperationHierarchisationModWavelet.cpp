/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/basis/modpoly/ModifiedPolyBasis.hpp"
#include "base/basis/modwavelet/operation/OperationHierarchisationModWavelet.hpp"



#include "base/exception/operation_exception.hpp"

namespace sg {
  namespace base {

    void OperationHierarchisationModWavelet::doHierarchisation(DataVector& node_values) {
      throw new operation_exception("This operation is not implemented, yet! Sorry ;-)");
    }

    void OperationHierarchisationModWavelet::doDehierarchisation(DataVector& alpha) {
      throw new operation_exception("This operation is not implemented, yet! Sorry ;-)");
    }

  }
}
