/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef ELLIPTICPDESOLVER_HPP
#define ELLIPTICPDESOLVER_HPP

#include "pde/application/PDESolver.hpp"

namespace sg {
  namespace pde {

    /**
     * This class extends the PDESolver with functions that are needed to
     * solve elliptic PDEs
     *
     * @version $HEAD$
     */
    class EllipticPDESolver : public PDESolver {
      protected:

      public:
        /**
         * Std-Constructor of the solver
         */
        EllipticPDESolver();

        /**
         * Std-Destructor of the solver
         */
        virtual ~EllipticPDESolver();

        /**
         * abstract method to solve an elliptic PDE. All solver of elliptic PDEs
         * have to implement this method.
         *
         * @param alpha the coefficients of the Sparse Gird's basis functions will be in this vector after solving
         * @param rhs the right hand side of the SLE
         * @param maxCGIterations the maximum of interation in the CG solver
         * @param epsilonCG the epsilon used in the CG
         * @param verbose enables verbose output during solving
         */
        virtual void solvePDE(sg::base::DataVector& alpha, sg::base::DataVector& rhs, size_t maxCGIterations, double epsilonCG, bool verbose = false) = 0;
    };

  }
}

#endif /* ELLIPTICPDESOLVER_HPP */
