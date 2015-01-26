/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef BICGSTABMPI_HPP
#define BICGSTABMPI_HPP

#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/base/operation/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <iostream>

namespace sg {
  namespace parallel {

    class BiCGStabMPI : public sg::solver::SLESolver {
      private:
        /**
         * Routine called by the MPI slaves, here just the execution of
         * of sub part of the SystemMatrix's mult-Routine is needed.
         *
         * @param SystemMatrix reference to an OperationMatrix Object that implements the matrix vector multiplication
         * @param alpha the sparse grid's coefficients which have to be determined
         */
        virtual void waitForTask(sg::base::OperationMatrix& SystemMatrix, sg::base::DataVector& alpha);

      public:
        /**
         * Std-Constructor
         */
        BiCGStabMPI(size_t imax, double epsilon);

        /**
         * Std-Destructor
         */
        virtual ~BiCGStabMPI();

        /**
         * max_threashold is ignored in this solver
         *
         * Reference:
         * http://www.iue.tuwien.ac.at/phd/heinreichsberger/node70.html
         * http://www.numerik.math.tu-graz.ac.at/kurse/lgs/SIMNET6.pdf
         * http://netlib.org
         */
        virtual void solve(sg::base::OperationMatrix& SystemMatrix, sg::base::DataVector& alpha, sg::base::DataVector& b, bool reuse = false, bool verbose = false, double max_threshold = -1.0);
    };

  }
}

#endif /*BICGSTABMPI_HPP */
