// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONIDENTITY_HPP
#define OPERATIONIDENTITY_HPP

#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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
        virtual ~OperationIdentity() override {}

        void mult(DataVector& alpha, DataVector& result) override {
          result = alpha;
        }
    };

  }
}
#endif /* OPERATIONIDENTITY_HPP */
