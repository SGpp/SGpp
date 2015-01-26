/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#ifndef OPERATIONIDENTITY_HPP
#define OPERATIONIDENTITY_HPP

#include <sgpp/base/operation/OperationMatrix.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>

namespace sg {
  namespace base {

    /**
     * Implementation of identity Operation for all kinds of grids
     */
    class OperationIdentity: public OperationMatrix {
      public:
        /**
         * Constructor of OperationIdentity
         */
        OperationIdentity() {
        }

        /**
         * Destructor
         */
        virtual ~OperationIdentity() {}

        void mult(DataVector& alpha, DataVector& result) {
          result = alpha;
        }
    };

  }
}
#endif /* OPERATIONIDENTITY_HPP */
