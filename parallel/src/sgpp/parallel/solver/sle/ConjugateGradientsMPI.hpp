/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef CONJUGATEGRADIENTSMPI_HPP
#define CONJUGATEGRADIENTSMPI_HPP

#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/base/operation/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    class ConjugateGradientsMPI : public SGPP::solver::SLESolver {
      private:
        /**
         * Routine called by the MPI slaves, here just the execution of
         * of sub part of the SystemMatrix's mult-Routine is needed.
         *
         * @param SystemMatrix reference to an SGPP::base::OperationMatrix Object that implements the matrix vector multiplication
         * @param alpha the sparse grid's coefficients which have to be determined
         */
        virtual void waitForTask(SGPP::base::OperationMatrix& SystemMatrix, SGPP::base::DataVector& alpha);

      public:
        ConjugateGradientsMPI(size_t imax = 0, double epsilon = 0.0);

        virtual ~ConjugateGradientsMPI();

        virtual void solve(SGPP::base::OperationMatrix& SystemMatrix, SGPP::base::DataVector& alpha, SGPP::base::DataVector& b, bool reuse = false, bool verbose = false, double max_threshold = -1.0);
    };

  }
}

#endif /* CONJUGATEGRADIENTSMPI_HPP */
