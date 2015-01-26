/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef OPERATIONQUADRATURE_HPP
#define OPERATIONQUADRATURE_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class provides the quadrature of a sparse grid function
     */
    class OperationQuadrature {
      public:
        /**
         * Constructor
         */
        OperationQuadrature() {}

        /**
         * Destructor
         */
        virtual ~OperationQuadrature() {}

        /**
         * Integrate the sparse grid function
         *
         * @param alpha the function's values in the nodal basis
         */
        virtual double doQuadrature(DataVector& alpha) = 0;

    };

  }
}

#endif /* OPERATIONQUADRATURE_HPP */
