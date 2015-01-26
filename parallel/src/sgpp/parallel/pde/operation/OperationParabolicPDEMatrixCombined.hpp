/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONPARABOLICPDEMATRIXCOMBINED_HPP
#define OPERATIONPARABOLICPDEMATRIXCOMBINED_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/OperationMatrix.hpp>

namespace sg {
  namespace parallel {

    /**
     * Abstract definition of a matrix operator interface used for
     * solving parabolic PDEs. It allows for applying both, the mass and stiffness
     * matrix in one mult-call.
     */
    class OperationParabolicPDEMatrixCombined : public sg::base::OperationMatrix {
      protected:
        /// storing the current timestep coefficient
        double TimestepCoeff;

      public:
        /**
         * Constructor
         */
        OperationParabolicPDEMatrixCombined() {}

        /**
         * Destructor
         */
        virtual ~OperationParabolicPDEMatrixCombined() {}

        /**
         * Sets the timestep coefficient
         *
         * @param newTimestepCoeff The new timestep coefficient for the chosen
         * numerical approximation scheme.
         */
        void setTimestepCoeff(double newTimestepCoeff) {
          this->TimestepCoeff = newTimestepCoeff;
        }

        /**
         * Gets the timestep coefficient
         *
         * @return newTimestepCoeff The new timestep coefficient for the chosen
         * numerical approximation scheme.
         */
        double getTimestepCoeff() {
          return this->TimestepCoeff;
        }

    };

  }
}

#endif /* OPERATIONPARABOLICPDEMATRIXCOMBINED_HPP */
