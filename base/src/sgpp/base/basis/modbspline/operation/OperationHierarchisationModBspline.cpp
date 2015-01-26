/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#include <sgpp/base/basis/modbspline/operation/OperationHierarchisationModBspline.hpp>

#include <sgpp/base/exception/operation_exception.hpp>

namespace sg {
  namespace base {

    void OperationHierarchisationModBspline::doHierarchisation(DataVector& node_values) {
      throw new operation_exception("This operation is not implemented, yet! Sorry ;-)");
    }

    void OperationHierarchisationModBspline::doDehierarchisation(DataVector& alpha) {
      throw new operation_exception("This operation is not implemented, yet! Sorry ;-)");
    }

  }
}
